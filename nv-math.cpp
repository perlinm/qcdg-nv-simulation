#include <iostream> // for standard output
using namespace std;

#define EIGEN_USE_MKL_ALL
#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "constants.h"
#include "qp-math.h"
#include "nv-math.h"

// ---------------------------------------------------------------------------------------
// Spin vectors and structs
// ---------------------------------------------------------------------------------------

// rotate into one basis from another
Matrix2cd rotate(const vector<Vector3d>& basis_end, const vector<Vector3d>& basis_start){
  assert(basis_start.size() == 3);
  assert(basis_end.size() == 3);
  assert(dot(basis_start.at(0).cross(basis_start.at(1)),basis_start.at(2)) > 0);
  assert(dot(basis_end.at(0).cross(basis_end.at(1)),basis_end.at(2)) > 0);

  // rotation matrix taking starting basis vectors to ending basis vectors
  const Matrix3d rotation = (basis_end.at(0)*basis_start.at(0).transpose() +
                             basis_end.at(1)*basis_start.at(1).transpose() +
                             basis_end.at(2)*basis_start.at(2).transpose());
  // rotation angle
  const double angle = acos((rotation.trace()-1.)/2.);

  // determine rotation axis
  const double axis_x = rotation(2,1)-rotation(1,2);
  const double axis_y = rotation(0,2)-rotation(2,0);
  const double axis_z = rotation(1,0)-rotation(0,1);
  const Vector3d axis = hat((Vector3d() << axis_x, axis_y, axis_z).finished());
  if(axis.squaredNorm() > 0){
    return rotate(angle, axis);
  } else{
    const EigenSolver<Matrix3d> solver(rotation);
    const Vector3cd e_vals = solver.eigenvalues();
    const Matrix3cd e_vecs = solver.eigenvectors();
    uint axis_index = 0;
    for(uint i = 1; i < e_vals.size(); i++){
      if(abs(e_vals(i)-1.) < abs(e_vals(axis_index)-1.)) axis_index = i;
    }
    return rotate(angle, e_vecs.col(axis_index).real());
  }
}

spin::spin(const Vector3d& pos, const double g, const mvec& S) :
  pos(pos), g(g), S(S)
{};

nv_system::nv_system(const int ms, const double static_Bz, const axy_harmonic k_DD,
                     const double scale_factor, const double integration_factor,
                     const bool no_nn) :
  e(Vector3d::Zero(), ge,
    mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat)),
  ms(ms), static_Bz(static_Bz), k_DD(k_DD),
  scale_factor(scale_factor), integration_factor(integration_factor), no_nn(no_nn)
{};

// ---------------------------------------------------------------------------------------
// Spin placement and clustering methods
// ---------------------------------------------------------------------------------------

// check whether a nucleus is addressable
bool can_address(const nv_system& nv, const uint target){
  const Vector3d r = nv.nuclei.at(target).pos - nv.e.pos;
  const Vector3i r_xy = xy_int_pos(r);
  const Vector3i r_z = z_int_pos(r);

  if(r_xy.squaredNorm() == 0 || r_z.squaredNorm() == 0){
    // target is on z axis or in x-y plane
    return false;
  }

  for(uint i = 0; i < nv.nuclei.size(); i++){
    if(i == target) continue;
    const Vector3d s = nv.nuclei.at(i).pos - nv.e.pos;
    const Vector3i s_xy = xy_int_pos(s);
    const Vector3i s_z = z_int_pos(s);
    if(s_z.squaredNorm() == r_z.squaredNorm() && (s_xy == r_xy || s_xy == -r_xy)){
      // target has a larmor pair with a parallel x-y component of the hyperfine field
      return false;
    }
  }

  return true;
}

// determine whether two spins are a larmor pair
bool is_larmor_pair(const nv_system& nv, const uint idx1, const uint idx2){
  const Vector3d r1 = nv.nuclei.at(idx1).pos - nv.e.pos;
  const Vector3d r2 = nv.nuclei.at(idx2).pos - nv.e.pos;
  return
    z_int_pos(r1).squaredNorm() == z_int_pos(r2).squaredNorm() &&
    xy_int_pos(r1).squaredNorm() == xy_int_pos(r2).squaredNorm();
}

// coupling strength between two spins; assumes strong magnetic field in zhat
double coupling_strength(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos - s1.pos;
  return abs(s1.g*s2.g/(8*pi*pow(r.norm()*a0/2,3)) *
             (1-3*dot(hat(r),zhat)*dot(hat(r),zhat)));
}

// group nuclei into clusters with intercoupling strengths >= coupling_strength
vector<vector<uint>> cluster_with_coupling(const nv_system& nv,
                                           const double min_coupling_strength,
                                           const bool cluster_by_larmor_frequency){

  vector<vector<uint>> clusters(0); // all clusters
  vector<uint> clustered(0); // indices of nuclei we have already clustered

  for(uint n_seed = 0; n_seed < nv.nuclei.size(); n_seed++){
    if(in_vector(n_seed,clustered)) continue;

    // initialize current cluster, and add the current nucleus
    vector<uint> cluster;
    cluster.push_back(n_seed);
    clustered.push_back(n_seed);

    // loop over all nuclei in cluster
    for(uint i = 0; i < cluster.size(); i++){
      const uint n_old = cluster.at(i);
      // loop over all nuclei, checking whether to add them to the cluster
      for(uint n_new = 0; n_new < nv.nuclei.size(); n_new++){
        if(in_vector(n_new,clustered)) continue;

        // if n_old and n_new are interacting or are a larmor pair,
        //   add n_new to this cluster
        if((coupling_strength(nv,n_old,n_new) > min_coupling_strength) ||
           (cluster_by_larmor_frequency && is_larmor_pair(nv,n_old,n_new))){
          cluster.push_back(n_new);
          clustered.push_back(n_new);
        }
      }
    }
    clusters.push_back(cluster);
  }
  return clusters;
}

// group together clusters sharing larmor pairs
vector<vector<uint>> group_clusters(const nv_system& nv){
  vector<vector<uint>> old_clusters = nv.clusters;
  vector<vector<uint>> new_clusters;

  while(old_clusters.size() > 0){

    new_clusters.push_back(old_clusters.at(0));
    old_clusters.erase(old_clusters.begin());

    bool grouped = false;
    vector<uint> new_cluster = new_clusters.back();
    for(uint i = 0; i < new_cluster.size(); i++){
      for(uint c = 0; c < old_clusters.size(); c++){
        const vector<uint> old_cluster = old_clusters.at(c);
        for(uint j = 0; j < old_clusters.at(c).size(); j++){
          if(is_larmor_pair(nv, new_cluster.at(i), old_cluster.at(j))){
            new_cluster.insert(new_cluster.end(), old_cluster.begin(), old_cluster.end());
            old_clusters.erase(old_clusters.begin()+c);
            c--;
            grouped = true;
            break;
          }
        }
      }
    }
    if(grouped){
      new_clusters.at(new_clusters.size()-1) = new_cluster;
    }
  }

  return new_clusters;
}

// get size of largest cluster
uint largest_cluster_size(const vector<vector<uint>>& clusters){
  uint largest_size = 0;
  for(uint i = 0; i < clusters.size(); i++){
    if(clusters.at(i).size() > largest_size) largest_size = clusters.at(i).size();
  }
  return largest_size;
}

// largest internuclear coupling
double largest_coupling(const nv_system& nv){
    double max = 0;
    for(uint i = 0; i < nv.nuclei.size(); i++){
      for(uint j = i+1; j < nv.nuclei.size(); j++){
        double c_ij = coupling_strength(nv,i,j);
        if(c_ij > max) max = c_ij;
      }
    }
    return max;
}

// cluster nuclei and set cluster_coupling for a given maximum cluster size
void cluster_nuclei(nv_system& nv, const uint max_cluster_size,
                    const bool cluster_by_larmor_frequency,
                    const double initial_cluster_coupling, const double cc_resolution){
  assert(cc_resolution > 0);
  uint cluster_size_target = max_cluster_size;

  // special cases for small clusters
  if(cluster_size_target >= nv.nuclei.size()){
    nv.cluster_coupling = largest_coupling(nv);
    nv.clusters.push_back({});
    for(uint n = 0; n < nv.nuclei.size(); n++){
      nv.clusters.at(0).push_back(n);
    }
    return;
  }
  const uint min_cluster_size_cap = smallest_possible_cluster_size(nv);
  if(cluster_size_target < min_cluster_size_cap){
    nv.cluster_coupling = largest_coupling(nv);
    for(uint n = 0; n < nv.nuclei.size(); n++){
      nv.clusters.push_back({n});
    }
    return;
  } else if(cluster_size_target == min_cluster_size_cap){
    cluster_size_target += 1;
  }

  // find largest coupling for which the largest cluster size is <= cluster_size_target
  nv.cluster_coupling = initial_cluster_coupling;
  double dcc = nv.cluster_coupling/8;

  bool coupling_too_small = true;
  bool last_coupling_too_small = true;
  bool crossed_correct_coupling = false;

  while(dcc >= cc_resolution){
    if(coupling_too_small == last_coupling_too_small){
      if(!crossed_correct_coupling) dcc *= 2;
    } else{
      dcc /= 2;
    }
    last_coupling_too_small = coupling_too_small;

    nv.cluster_coupling += coupling_too_small ? dcc : -dcc;
    nv.clusters = cluster_with_coupling(nv, nv.cluster_coupling,
                                        cluster_by_larmor_frequency);

    coupling_too_small = (largest_cluster_size(nv.clusters) > cluster_size_target);
    if(!crossed_correct_coupling && (coupling_too_small != last_coupling_too_small)){
      crossed_correct_coupling = true;
    }
  }

  while(largest_cluster_size(nv.clusters) > cluster_size_target){
    nv.cluster_coupling += dcc;
    nv.clusters = cluster_with_coupling(nv, nv.cluster_coupling,
                                        cluster_by_larmor_frequency);
  }
}

uint get_cluster_containing_target(const nv_system& nv, const uint index){
  assert(index < nv.nuclei.size());
  for(uint c = 0; c < nv.clusters.size(); c++){
    for(uint s = 0; s < nv.clusters.at(c).size(); s++){
      if(nv.clusters.at(c).at(s) == index){
        return c;
      }
    }
  }
  return 0;
}

uint get_index_in_cluster(const uint index, const vector<uint> cluster){
  assert(in_vector(index,cluster));
  for(uint s = 0; s < cluster.size(); s++){
    if(cluster.at(s) == index){
      return s;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------------------
// AXY scanning methods
// ---------------------------------------------------------------------------------------

// component of hyperfine field perpendicular to the larmor axis
Vector3d hyperfine_perp(const nv_system&nv, const spin& s){
  const Vector3d w_eff = effective_larmor(nv,s);
  const Vector3d A = hyperfine(nv,s);
  return A - dot(A,hat(w_eff))*hat(w_eff);
}

// minimum difference in larmor frequencies between target nucleus and other nuclei
//   i.e. min{ |w_s - w_{index}| for all s != index }
double larmor_resolution(const nv_system& nv, const uint index){
  const double target_larmor = effective_larmor(nv,index).norm();
  double dw_min = target_larmor; // maximum allowable larmor resolution
  for(uint s = 0; s < nv.nuclei.size(); s++){
    if(is_larmor_pair(nv,s,index)) continue; // find dw_min for distinct frequencies only
    const double dw = abs(target_larmor - effective_larmor(nv,s).norm());
    if(dw < dw_min) dw_min = dw;
  }
  return dw_min;
}

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const double f, const axy_harmonic k){
  assert(abs(f) <= axy_f_max(k));
  const double fp = f*pi;

  // compute first two pulse times
  double x1,x2;
  if(k == 1){
    const double w1 = 4 - fp;
    const double w2 = w1 * (960 - 144*fp - 12*fp*fp + fp*fp*fp);

    x1 = 1/(2*pi) * atan2( (3*fp-12)*w1 + sqrt(3*w2),
                           sqrt(6)*sqrt(w2 - 96*fp*w1 + w1*w1*sqrt(3*w2)));
    x2 = 1/(2*pi) * atan2(-(3*fp-12)*w1 + sqrt(3*w2),
                          sqrt(6)*sqrt(w2 - 96*fp*w1 - w1*w1*sqrt(3*w2)));
  } else{ // if k == 3
    const double q1 = 4/(sqrt(5+fp)-1);
    const double q2 = 4/(sqrt(5+fp)+1);

    x1 = 1./4 - 1/(2*pi)*atan(sqrt(q1*q1-1));
    x2 = 1./4 - 1/(2*pi)*atan(sqrt(q2*q2-1));
  }

  // construct vector of all pulse times (normalized to one AXY period)
  vector<double> pulse_times;
  pulse_times.push_back(0);
  pulse_times.push_back(x1);
  pulse_times.push_back(x2);
  pulse_times.push_back(0.25);
  pulse_times.push_back(0.5 - x2);
  pulse_times.push_back(0.5 - x1);
  pulse_times.push_back(0.5 + x1);
  pulse_times.push_back(0.5 + x2);
  pulse_times.push_back(0.75);
  pulse_times.push_back(1. - x2);
  pulse_times.push_back(1. - x1);
  pulse_times.push_back(1);
  return pulse_times;
}

vector<double> advanced_pulse_times(const vector<double> pulse_times,
                                    const double advance){
  const double normed_advance = mod(advance, 1);
  if(normed_advance == 0) return pulse_times;

  // number of pulses
  const uint N = pulse_times.size()-2;

  // advanced pulse_times
  vector<double> advanced_pulse_times;
  advanced_pulse_times.push_back(0);
  for(uint p = 0; p < 2*N; p++){
    if(p/N + pulse_times.at(p%N+1) - normed_advance >= 0){
      advanced_pulse_times.push_back(p/N + pulse_times.at(p%N+1) - normed_advance);
    }
    if(advanced_pulse_times.size()-1 == N) break;
  }
  advanced_pulse_times.push_back(1);
  return advanced_pulse_times;
}

// evaluate F(x) (i.e. sign in front of sigma_z^{NV}) for given AXY pulses
int F_AXY(const double x, const vector<double> pulses){
  const double normed_x = mod(x, 1);
  uint pulse_count = 0;
  for(uint i = 1; i < pulses.size()-1; i++){
    if(pulses.at(i) < normed_x) pulse_count++;
    else break;
  }
  return (pulse_count % 2 == 0 ? 1 : -1);
}

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos - s1.pos;
  return s1.g*s2.g/(4*pi*pow(r.norm()*a0/2,3))
    * (dot(s1.S,s2.S) - 3*tp(dot(s1.S,hat(r)), dot(s2.S,hat(r))));
}

// spin coupling Hamiltonian for the entire spin system
MatrixXcd H_int(const nv_system& nv, const uint cluster_index){
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction between NV center and spin s
    const spin n_s = nv.nuclei.at(cluster.at(s));
    H += act(H_ss(nv.e,n_s), {0,s+1}, spins);
    if(!nv.no_nn){
      // interaction betwen spin r and spin s
      for(uint r = 0; r < s; r++){
        const spin n_r = nv.nuclei.at(cluster.at(r));
        H += act(H_ss(n_s,n_r), {s+1,r+1}, spins);
      }
    }
  }
  return H;
}

// nuclear Zeeman Hamiltonian
MatrixXcd H_nZ(const nv_system& nv, const uint cluster_index, const Vector3d& B){
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  const int spins = cluster.size()+1;
  // zero-field splitting and interaction of NV center with magnetic field
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction of spin s with magnetic field
    H -= act(nv.nuclei.at(cluster.at(s)).g*dot(B,nv.nuclei.at(cluster.at(s)).S),
             {s+1}, spins);
  }
  return H;
}

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time){
  const double w_DD = w_scan/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulse_times = axy_pulse_times(f_DD,nv.k_DD); // AXY pulse times

  double coherence = 1;
  for(uint cluster = 0; cluster < nv.clusters.size(); cluster++){
    const uint cluster_size = nv.clusters.at(cluster).size();

    // projections onto |ms> and |0> NV states
    const MatrixXcd proj_m = act(up*up.adjoint(),{0},cluster_size+1); // |ms><ms|
    const MatrixXcd proj_0 = act(dn*dn.adjoint(),{0},cluster_size+1); // |0><0|

    // full system Hamiltonian
    const MatrixXcd H = H_sys(nv,cluster);

    // spin cluster Hamiltonians for single NV states
    const MatrixXcd H_m = ptrace(H*proj_m, {0}); // <ms|H|ms>
    const MatrixXcd H_0 = ptrace(H*proj_0, {0}); // <0|H|0>

    // propagators for sections of the AXY sequence
    const MatrixXcd U1_m = exp(-j*H_m*t_DD*(pulse_times.at(1)-pulse_times.at(0)));
    const MatrixXcd U2_m = exp(-j*H_0*t_DD*(pulse_times.at(2)-pulse_times.at(1)));
    const MatrixXcd U3_m = exp(-j*H_m*t_DD*(pulse_times.at(3)-pulse_times.at(2)));

    const MatrixXcd U1_0 = exp(-j*H_0*t_DD*(pulse_times.at(1)-pulse_times.at(0)));
    const MatrixXcd U2_0 = exp(-j*H_m*t_DD*(pulse_times.at(2)-pulse_times.at(1)));
    const MatrixXcd U3_0 = exp(-j*H_0*t_DD*(pulse_times.at(3)-pulse_times.at(2)));

    // single AXY sequence propagators
    MatrixXcd U_m = U1_m*U2_m*U3_m*U3_0*U2_0*U1_0 * U1_0*U2_0*U3_0*U3_m*U2_m*U1_m;
    MatrixXcd U_0 = U1_0*U2_0*U3_0*U3_m*U2_m*U1_m * U1_m*U2_m*U3_m*U3_0*U2_0*U1_0;

    // propagators for entire scan
    U_m = pow(U_m, int(scan_time/t_DD));
    U_0 = pow(U_0, int(scan_time/t_DD));

    // normalize propagators
    U_m /= sqrt(real(trace(U_m.adjoint()*U_m)/double(U_m.rows())));
    U_0 /= sqrt(real(trace(U_0.adjoint()*U_0)/double(U_0.rows())));

    // update coherence
    coherence *= real(trace(U_0.adjoint()*U_m)) / pow(2,cluster_size);
  }
  return coherence;
}

// ---------------------------------------------------------------------------------------
// Control fields and simulation
// ---------------------------------------------------------------------------------------

// return control field for decoupling a single nucleus from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd){
  const spin s = nv.nuclei.at(index);
  const Vector3d w_j = effective_larmor(nv,index);
  const double w_rfd = w_j.norm()/(1-sin(theta_rfd)/(2*sqrt(2)*nv.scale_factor));
  const double V_rfd = w_rfd/(s.g*nv.scale_factor);
  const Vector3d n_rfd = rotate(hat(hyperfine_perp(nv, index)), theta_rfd, hat(w_j));
  return control_fields(V_rfd*n_rfd, w_rfd, phi_rfd);
}

// perform given rotation on the NV center
MatrixXcd rotate_NV(const nv_system& nv, const Vector3d& rotation, const uint spins){
  if(rotation.squaredNorm() > 0){
    return act( exp(-j*dot(rotation,nv.e.S)), {0}, spins);
  } else{
    return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  }
}

// compute and perform rotation of NV center necessary to generate U_NV
MatrixXcd act_NV(const nv_system& nv, const Matrix2cd& U_NV, const uint spins){
  const Vector4cd H_NV_vec = U_decompose(j*log(U_NV));
  const Vector3d nv_rotation = (xhat*real(H_NV_vec(1))*sqrt(2) +
                                yhat*real(H_NV_vec(2))*nv.ms*sqrt(2) +
                                zhat*real(H_NV_vec(3))*nv.ms*2);
  return rotate_NV(nv,nv_rotation,spins);
}

// simulate propagator with static control fields
MatrixXcd simulate_AXY8(const nv_system& nv, const uint cluster,
                        const double w_DD, const double f_DD, const axy_harmonic k_DD,
                        const double simulation_time, const double advance,
                        const Vector3d B_ctl){
  const uint spins = nv.clusters.at(cluster).size()+1;

  // AXY sequence parameters
  const double t_DD = 2*pi/w_DD;
  const vector<double> pulses = axy_pulse_times(f_DD,k_DD);
  const vector<double> advanced_pulses = advanced_pulse_times(pulses, advance/t_DD);

  const MatrixXcd H = H_sys(nv,cluster) + H_ctl(nv,cluster,B_ctl); // full Hamiltonian
  const MatrixXcd X = act_NV(nv, sx, spins); // NV center spin flip (pi-)pulse
  MatrixXcd U = MatrixXcd::Identity(H.rows(),H.cols()); // system propagator
  Matrix2cd U_NV = I2; // NV-only propagator

  // if we need to start with a flipped NV center, flip it
  if(F_AXY(advance/t_DD, pulses) == -1){
    U = (X * U).eval();
    U_NV = (sx * U_NV).eval();
  }

  // propagator for whole AXY sequences
  if(simulation_time >= t_DD){
    MatrixXcd U_AXY = MatrixXcd::Identity(H.rows(),H.cols());
    Matrix2cd U_NV_AXY = I2;
    for(uint i = 1; i < advanced_pulses.size(); i++){
      const double t = advanced_pulses.at(i-1)*t_DD;
      const double dt = advanced_pulses.at(i)*t_DD - t;
      U_AXY = (X * exp(-j*dt*H) * U_AXY).eval();
      U_NV_AXY = (sx * exp(-j*dt*H_NV(nv,B_ctl)) * U_NV_AXY).eval();
    }
    // "undo" the last pulse at t = t_DD
    U_AXY = (X * U_AXY).eval();
    U_NV_AXY = (sx * U_NV_AXY).eval();

    U = (pow(U_AXY, int(simulation_time/t_DD)) * U).eval();
    U_NV = (pow(U_NV_AXY, int(simulation_time/t_DD)) * U_NV).eval();
  }

  // propagator for AXY sequence remainder
  const double remaining_time = simulation_time - int(simulation_time/t_DD)*t_DD;
  for(uint i = 1; i < advanced_pulses.size(); i++){
    const double t = advanced_pulses.at(i-1)*t_DD;
    const double dt = advanced_pulses.at(i)*t_DD - t;
    if(t + dt < remaining_time){
      U = (X * exp(-j*dt*H) * U).eval();
      U_NV = (sx * exp(-j*dt*H_NV(nv,B_ctl)) * U_NV).eval();
    } else{
      const double dtf = remaining_time - t;
      U = (exp(-j*dtf*H) * U).eval();
      U_NV = (exp(-j*dtf*H_NV(nv,B_ctl)) * U_NV).eval();
      break;
    }
  }
  // rotate into the frame of the NV center and normalize the propagator
  U = (act_NV(nv,U_NV.adjoint(),spins) * U).eval();
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));
  return U;
}

// simulate propagator with dynamic control fields
MatrixXcd simulate_AXY8(const nv_system& nv, const uint cluster,
                        const double w_DD, const double f_DD, const axy_harmonic k_DD,
                        const control_fields& controls, const double simulation_time,
                        const double advance){
  if(controls.all_fields_static()){
    return simulate_AXY8(nv, cluster, w_DD, f_DD, k_DD,
                         simulation_time, advance, controls.B(0));
  }
  const uint spins = nv.clusters.at(cluster).size()+1;
  if(simulation_time == 0) return MatrixXcd::Identity(pow(2,spins),pow(2,spins));

  // AXY sequence parameters
  const double t_DD = 2*pi/w_DD;
  const vector<double> pulses = axy_pulse_times(f_DD,k_DD);

  // largest frequency scale of simulation
  const double frequency_scale = [&]() -> double {
    double largest_control_freq = w_DD;
    Vector3d B_cap = abs(nv.static_Bz)*zhat;
    for(uint c = 0; c < controls.num(); c++){
      largest_control_freq = max(largest_control_freq,controls.freqs.at(c));
      B_cap += (abs(dot(controls.Bs.at(c),xhat))*xhat +
                abs(dot(controls.Bs.at(c),yhat))*yhat +
                abs(dot(controls.Bs.at(c),zhat))*zhat);
    }
    double largest_g = 0;
    for(uint n = 0; n < nv.clusters.at(cluster).size(); n++){
      largest_g = max(largest_g,abs(nv.nuclei.at(nv.clusters.at(cluster).at(n)).g));
    }
    return max(largest_control_freq,largest_g*B_cap.norm());
  }();

  // integration step number and size
  const uint integration_steps = simulation_time*frequency_scale*nv.integration_factor;
  const double dt = simulation_time/integration_steps;

  const MatrixXcd H_0 = H_sys(nv,cluster); // full system Hamiltonian
  const MatrixXcd X = act_NV(nv, sx, spins); // NV center spin flip (pi-)pulse
  MatrixXcd U = MatrixXcd::Identity(H_0.rows(),H_0.cols()); // system propagator
  Matrix2cd U_NV = I2; // NV-only propagator

  // if we need to start with a flipped NV center, flip it
  if(F_AXY(advance/t_DD, pulses) == -1){
    U = (X * U).eval();
    U_NV = (sx * U_NV).eval();
  }

  for(uint t_i = 0; t_i < integration_steps; t_i++){
    const double t = t_i*dt + advance; // time

    // determine whether to apply an NV pi-pulse
    uint pulse = 0;
    double t_AXY = mod(t, t_DD); // time into this AXY sequence
    for(uint p = 1; p < pulses.size()-1; p++){
      if(pulses.at(p)*t_DD < t_AXY + dt){
        if(pulses.at(p)*t_DD >= t_AXY){
          pulse = p;
          break;
        }
      } else break;
    }

    // update propagator
    if(!pulse){
      const Vector3d B = controls.B(t+dt/2);
      const MatrixXcd H = H_0 + H_ctl(nv, cluster, B);

      U = (exp(-j*dt*H) * U).eval();
      U_NV = (exp(-j*dt*H_NV(nv,B)) * U_NV).eval();

    } else{ // if(pulse)
      const double t_AXY_initial = t_AXY;
      const double dt_full = dt;

      bool overflow = false;
      do{
        const double dt = (pulses.at(pulse)+overflow)*t_DD - t_AXY; // time before pulse
        const Vector3d B = controls.B(t+dt/2);
        const MatrixXcd H = H_0 + H_ctl(nv, cluster, B);

        U = (X * exp(-j*dt*H) * U).eval();
        U_NV = (sx * exp(-j*dt*H_NV(nv,B)) * U_NV).eval();

        t_AXY = (pulses.at(pulse)+overflow)*t_DD;
        pulse++;
        if(pulse == pulses.size()-1){
          pulse = 1;
          overflow = true;
        }
      } while((pulses.at(pulse)+overflow)*t_DD - t_AXY_initial < dt_full);

      const double dt = t_AXY_initial + dt_full - t_AXY; // time after last pulse
      const Vector3d B = controls.B(t+dt/2);
      const MatrixXcd H = H_0 + H_ctl(nv, cluster, B);

      U = (exp(-j*dt*H) * U).eval();
      U_NV = (exp(-j*dt*H_NV(nv,B)) * U_NV).eval();
    }
  }
  // rotate into the frame of the NV center and normalize the propagator
  U = (act_NV(nv,U_NV.adjoint(),spins) * U).eval();
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));
  return U;
}
