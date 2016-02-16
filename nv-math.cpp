#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "constants.h"
#include "qp-math.h"
#include "nv-math.h"

//--------------------------------------------------------------------------------------------
// Spin vectors and structs
//--------------------------------------------------------------------------------------------

// perform spin-1/2 rotation about arbitrary axis
Matrix2cd rotate(const Vector3d axis, const double phi){
  if(axis.squaredNorm() > 0) return cos(phi/2.)*I2 - j*sin(phi/2.)*dot(s_vec,hat(axis));
  else return I2;
}

// rotate into one basis from another
Matrix2cd rotate(const vector<Vector3d> basis_end, const vector<Vector3d> basis_start){
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
    return rotate(axis,angle);
  } else{
    const EigenSolver<Matrix3d> solver(rotation);
    const Vector3cd e_vals = solver.eigenvalues();
    const Matrix3cd e_vecs = solver.eigenvectors();
    uint axis_index = 0;
    for(uint i = 1; i < e_vals.size(); i++){
      if(abs(e_vals(i)-1.) < abs(e_vals(axis_index)-1.)) axis_index = i;
    }
    return rotate(e_vecs.col(axis_index).real(),angle);
  }
}

spin::spin(const Vector3d pos, const double g, const mvec S) :
  pos(pos), g(g), S(S)
{};

nv_system::nv_system(const int ms, const double static_Bz,
                     const double scale_factor, const uint integration_factor) :
  n(ao, 0., s_vec/2),
  e(spin(Vector3d::Zero(), ge,
         mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat))),
  ms(ms), static_Bz(static_Bz),
  scale_factor(scale_factor), integration_factor(integration_factor)
{};

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// determine whether two spins are a larmor pair
bool is_larmor_pair(const nv_system& nv, const uint idx1, const uint idx2){
  const Vector3d r1 = nv.nuclei.at(idx1).pos - nv.e.pos;
  const Vector3d r2 = nv.nuclei.at(idx2).pos - nv.e.pos;

  const int par_1 = round(16*abs(dot(r1,ao)));
  const int par_2 = round(16*abs(dot(r2,ao)));

  const int perp_1 = round(12*(r1-dot(r1,zhat)*zhat).squaredNorm());
  const int perp_2 = round(12*(r2-dot(r2,zhat)*zhat).squaredNorm());

  return par_1 == par_2 && perp_1 == perp_2;
}

// coupling strength between two spins; assumes strong magnetic field in zhat
inline double coupling_strength(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos - s1.pos;
  return abs(s1.g*s2.g/(8*pi*pow(r.norm()*a0,3)) * (1 - 3*dot(hat(r),zhat)*dot(hat(r),zhat)));
}

// group nuclei into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<uint>> cluster_nuclei(const vector<spin>& nuclei,
                                    const double min_coupling_strength){

  vector<vector<uint>> clusters; // all clusters
  vector<uint> clustered; // indices of nuclei we have already clustered

  // loop over indices of all nuclei
  for(uint i = 0; i < nuclei.size(); i++){
    if(!in_vector(i,clustered)){

      // initialize current cluster
      vector<uint> cluster;

      cluster.push_back(i);
      clustered.push_back(i);

      // loop over all indices of cluster
      for(uint ci = 0; ci < cluster.size(); ci++) {

        // loop over all nuclei indices greater than cluster(ci)
        for(uint k = cluster.at(ci)+1; k < nuclei.size(); k++){

          // if cluster(ci) and nuclei(k) are interacting, add k to this cluster
          if(!in_vector(k,clustered) &&
             coupling_strength(nuclei.at(cluster.at(ci)),
                               nuclei.at(k)) >= min_coupling_strength){
            cluster.push_back(k);
            clustered.push_back(k);
          }
        }
      }
      clusters.push_back(cluster);
    }
  }
  return clusters;
}

// group together clusters close nuclei have similar larmor frequencies
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

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(const vector<spin>& nuclei, const double initial_cluster_coupling,
                            const uint cluster_size_target,  const double dcc_cutoff){
  assert(dcc_cutoff > 0);

  // special case for "clusters" of 1 nucleus
  if(cluster_size_target == 1){
    double max_coupling = 0;
    for(uint i = 0; i < nuclei.size(); i++){
      for(uint j = i+1; j < nuclei.size(); j++){
        double c_ij = coupling_strength(nuclei.at(i), nuclei.at(j));
        if(c_ij > max_coupling) max_coupling = c_ij;
      }
    }
    return max_coupling + dcc_cutoff;
  }

  double cluster_coupling = initial_cluster_coupling;
  double dcc = cluster_coupling/4;

  vector<vector<uint>> clusters = cluster_nuclei(nuclei, cluster_coupling);
  bool coupling_too_small = largest_cluster_size(clusters) >= cluster_size_target;
  bool last_coupling_too_small;
  bool crossed_correct_coupling = false;

  // determine coupling for which the largest cluster size just barely >= max_cluster_size
  while(dcc >= dcc_cutoff || !coupling_too_small){
    last_coupling_too_small = coupling_too_small;

    cluster_coupling += coupling_too_small ? dcc : -dcc;
    clusters = cluster_nuclei(nuclei,cluster_coupling);
    coupling_too_small = largest_cluster_size(clusters) >= cluster_size_target;

    if(coupling_too_small != last_coupling_too_small) crossed_correct_coupling = true;

    if(coupling_too_small == last_coupling_too_small){
      if(!crossed_correct_coupling) dcc *= 2;
    } else{
      dcc /= 2;
    }
  }
  return cluster_coupling;
}

uint get_cluster_containing_index(const nv_system& nv, const uint index){
  assert(index < nv.nuclei.size());
  for(uint c = 0; c < nv.clusters.size(); c++){
    for(uint s = 0; s < nv.clusters.at(c).size(); s++){
      if(nv.clusters.at(c).at(s) == index){
        return c;
      }
    }
  }
}

uint get_index_in_cluster(const uint index, const vector<uint> cluster){
  assert(in_vector(index,cluster));
  for(uint s = 0; s < cluster.size(); s++){
    if(cluster.at(s) == index){
      return s;
    }
  }
}

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// component of hyperfine field perpendicular to the larmor axis
Vector3d A_perp(const nv_system&nv, const spin& s){
  const Vector3d larmor_eff = effective_larmor(nv,s);
  const Vector3d hyperfine = A(nv,s);
  return hyperfine - dot(hyperfine,hat(larmor_eff))*hat(larmor_eff);
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

vector<double> advanced_pulse_times(const vector<double> pulse_times, const double advance){
  const double normed_advance = advance - floor(advance);
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
    if(advanced_pulse_times.size() == N+1) break;
  }
  advanced_pulse_times.push_back(1);
  return advanced_pulse_times;
}

// evaluate F(x) (i.e. sign in front of sigma_z^{NV}) for given AXY pulses
int F_AXY(const double x, const vector<double> pulses){
  uint pulse_count = 0;
  for(uint i = 1; i < pulses.size()-1; i++){
    if(pulses.at(i) < x - floor(x)) pulse_count++;
    else break;
  }
  return (pulse_count%2 == 0 ? 1 : -1);
}

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos-s1.pos;
  return s1.g*s2.g/(4*pi*pow(r.norm()*a0,3))
    * (dot(s1.S,s2.S) - 3*tp(dot(s1.S,hat(r)), dot(s2.S,hat(r))));
}

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const nv_system& nv, const uint cluster_index){
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction between NV center and spin s
    H += act(H_ss(nv.e, nv.nuclei.at(cluster.at(s))), {0,s+1}, spins);
    for(uint r = 0; r < s; r++){
      // interaction betwen spin r and spin s
      H += act(H_ss(nv.nuclei.at(cluster.at(r)), nv.nuclei.at(cluster.at(s))),
               {r+1,s+1}, spins);
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

// Zeeman Hamiltonian for NV+cluster system
inline MatrixXcd H_Z(const nv_system& nv, const uint cluster, const Vector3d& B){
  return H_nZ(nv,cluster,B) + act(H_NV_GS(nv,B), {0}, nv.clusters.at(cluster).size()+1);
}

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const axy_harmonic k_DD, const double scan_time){
  const double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulse_times = axy_pulse_times(f_DD,k_DD); // AXY pulse times

  double coherence = 1;
  for(uint cluster = 0; cluster < nv.clusters.size(); cluster++){
    const uint cluster_size = nv.clusters.at(cluster).size();

    // projections onto |ms> and |0> NV states
    const MatrixXcd proj_m = act(up*up.adjoint(),{0},cluster_size+1); // |ms><ms|
    const MatrixXcd proj_0 = act(dn*dn.adjoint(),{0},cluster_size+1); // |0><0|

    // construct full Hamiltonian
    const MatrixXcd H = H_int(nv,cluster) + H_Z(nv,cluster,nv.static_Bz*zhat);

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

//--------------------------------------------------------------------------------------------
// Control fields and simulation
//--------------------------------------------------------------------------------------------

// return control field for decoupling a single nucleus from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd){
  const spin s = nv.nuclei.at(index);
  const Vector3d w_j = effective_larmor(nv,index);
  const double w_rfd = w_j.norm()/(1-sin(theta_rfd)/(2*sqrt(2)*nv.scale_factor));
  const double V_rfd = w_rfd/(s.g*nv.scale_factor);
  const Vector3d n_rfd = cos(theta_rfd)*hat(w_j) + sin(theta_rfd)*hat(w_j.cross(zhat));
  return control_fields(V_rfd*n_rfd, w_rfd, phi_rfd);
}

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const double f_DD, const axy_harmonic k_DD,
                              const double simulation_time, const double advance){
  const uint spins = nv.clusters.at(cluster).size()+1;

  // AXY sequence parameters
  const double t_DD = 2*pi/w_DD;
  const double normed_advance = advance/t_DD - floor(advance/t_DD);
  const vector<double> pulses = axy_pulse_times(f_DD,k_DD);
  const vector<double> advanced_pulses = advanced_pulse_times(pulses, normed_advance);

  // NV+cluster Hamiltonian
  const MatrixXcd H = H_int(nv,cluster) + H_Z(nv,cluster,nv.static_Bz*zhat);

  // NV center spin flip (pi-)pulse
  const MatrixXcd X = act(sx, {0}, spins);

  // initial propagator; determine whether to start with a flipped NV center (i.e. F(t) = -1)
  MatrixXcd U = MatrixXcd::Identity(H.rows(),H.cols());
  uint pulse_count = 0;
  for(uint i = 1; i < pulses.size()-1; i++){
    if(pulses.at(i) < normed_advance) pulse_count++;
    else break;
  }
  if(pulse_count%2 != 0) U = (X*U).eval();

  // propagator for whole AXY sequences
  if(simulation_time >= t_DD){
    MatrixXcd U_AXY = MatrixXcd::Identity(H.rows(),H.cols());
    for(uint i = 1; i < advanced_pulses.size(); i++){
      U_AXY = (X *
               U_NV_GS(nv, advanced_pulses.at(i)*t_DD, spins).adjoint() *
               exp(-j*H*(advanced_pulses.at(i)-advanced_pulses.at(i-1))*t_DD) *
               U_NV_GS(nv, advanced_pulses.at(i-1)*t_DD, spins) *
               U_AXY).eval();
    }
    U_AXY = (X * U_AXY).eval(); // "undo" the last pulse at t = t_DD
    U = (pow(U_AXY, int(simulation_time/t_DD)) * U).eval();
  }

  // propagator for AXY sequence remainder
  const double remaining_time = simulation_time - int(simulation_time/t_DD)*t_DD;
  for(uint i = 1; i < advanced_pulses.size(); i++){
    if(advanced_pulses.at(i)*t_DD < remaining_time){
      U = (X *
           U_NV_GS(nv, advanced_pulses.at(i)*t_DD, spins).adjoint() *
           exp(-j*H*t_DD*(advanced_pulses.at(i)-advanced_pulses.at(i-1))) *
           U_NV_GS(nv, advanced_pulses.at(i-1)*t_DD, spins) *
           U).eval();
      pulse_count++;
    } else{
      U = (U_NV_GS(nv, remaining_time, spins).adjoint() *
           exp(-j*H*(remaining_time-advanced_pulses.at(i-1)*t_DD)) *
           U_NV_GS(nv, advanced_pulses.at(i-1)*t_DD, spins) *
           U).eval();
      break;
    }
  }
  // if we ended with a flipped NV center (i.e. F(t) = -1), flip it back
  if(pulse_count%2 != 0) U = (X*U).eval();

  // normalize the propagator
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));

  return U;
}

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const double f_DD, const axy_harmonic k_DD,
                              const double simulation_time, const control_fields& controls,
                              const double advance){
  const uint spins = nv.clusters.at(cluster).size()+1;
  const double end_time = simulation_time + advance;

  // AXY sequence parameters
  const double t_DD = 2*pi/w_DD;
  const double normed_advance = advance/t_DD - floor(advance/t_DD);
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

  // static NV+cluster Hamiltonian
  const MatrixXcd H_static = H_int(nv,cluster);

  // NV center spin flip (pi-)pulse
  const MatrixXcd X = act(sx, {0}, spins);

  // initial propagator; determine whether to start with a flipped NV center (i.e. F(t) = -1)
  MatrixXcd U = MatrixXcd::Identity(H_static.rows(),H_static.cols());
  uint pulse_count = 0;
  for(uint i = 1; i < pulses.size()-1; i++){
    if(pulses.at(i) < normed_advance) pulse_count++;
    else break;
  }
  if(pulse_count%2 != 0) U = (X*U).eval();

  for(uint t_i = 0; t_i < integration_steps; t_i++){
    const double t = t_i*dt+advance; // time

    // determine whether to apply an NV pi-pulse
    uint pulse = 0;
    for(uint p = 1; p < pulses.size()-1; p++){
      if(pulses.at(p)*t_DD > t && pulses.at(p)*t_DD < t + dt){
        pulse = p;
        break;
      }
    }

    // update propagator
    if(!pulse){
      const Vector3d B = nv.static_Bz*zhat + controls.B(t+dt/2);
      const MatrixXcd H = H_static + H_nZ(nv,cluster,B);

      U = (exp(-j*dt*H) * U).eval();

    } else{ // if(pulse);
      const double t_AXY = t - int(t/t_DD)*t_DD; // time into this AXY sequence
      const double dt1 = pulses.at(pulse)*t_DD - t_AXY; // time before the pulse
      const double dt2 = dt - dt1; // time after the pulse

      const Vector3d B1 = nv.static_Bz*zhat + controls.B(t+dt1/2);
      const Vector3d B2 = nv.static_Bz*zhat + controls.B(t+dt1+dt2/2);
      const MatrixXcd H1 = H_static + H_nZ(nv,cluster,B1);
      const MatrixXcd H2 = H_static + H_nZ(nv,cluster,B2);

      U = (exp(-j*dt2*H2) * X * exp(-j*dt1*H1) * U).eval();
      pulse_count++;
    }
  }

  // if we ended with a flipped NV center (i.e. F(t) = -1), flip it back
  if(pulse_count%2 != 0) U = (X*U).eval();

  // normalize propagator
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));

  return U;
}
