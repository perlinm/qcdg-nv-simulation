#include <iostream> // for standard output
using namespace std;

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

nv_system::nv_system(const int ms, const double static_Bz, const double scale_factor) :
  n(ao, 0., s_vec/2),
  e(spin(Vector3d::Zero(), ge,
         mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat))),
  ms(ms), static_Bz(static_Bz), scale_factor(scale_factor)
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

uint get_index_in_cluster(const nv_system& nv, const uint index){
  assert(index < nv.nuclei.size());
  for(uint c = 0; c < nv.clusters.size(); c++){
    for(uint s = 0; s < nv.clusters.at(c).size(); s++){
      if(nv.clusters.at(c).at(s) == index){
        return s;
      }
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
  double dw_min = DBL_MAX;
  for(uint s = 0; s < nv.nuclei.size(); s++){
    if(s == index) continue;
    const double dw = abs(target_larmor - effective_larmor(nv,s).norm());
    if(dw < dw_min) dw_min = dw;
  }
  return dw_min;
}

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const uint k, const double f){

  assert((k == 1) || (k == 3));

  // compute first two pulse times
  const double fp = f*pi;
  double x1,x2;
  if(k == 1){
    assert(abs(fp) < 8*cos(pi/9)-4);
    const double w1 = 4 - fp;
    const double w2 = w1 * (960 - 144*fp - 12*fp*fp + fp*fp*fp);

    x1 = 1/(2*pi) * atan2( (3*fp-12)*w1 + sqrt(3*w2),
                           sqrt(6)*sqrt(w2 - 96*fp*w1 + w1*w1*sqrt(3*w2)));
    x2 = 1/(2*pi) * atan2(-(3*fp-12)*w1 + sqrt(3*w2),
                          sqrt(6)*sqrt(w2 - 96*fp*w1 - w1*w1*sqrt(3*w2)));
  } else{ // if k == 3
    assert(abs(fp) < 4);
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

vector<double> delayed_pulse_times(const vector<double> pulse_times, const double delay){
  // number of pulses
  const uint N = pulse_times.size()-2;

  // delayed pulse_times
  vector<double> delayed_pulse_times;
  delayed_pulse_times.push_back(0);
  for(uint p = 0; p < 2*N; p++){
    if(p/N + pulse_times.at(p%N+1) - delay >= 0){
      delayed_pulse_times.push_back(p/N + pulse_times.at(p%N+1) - delay);
    }
    if(delayed_pulse_times.size() == N+1) break;
  }
  delayed_pulse_times.push_back(1);
  return delayed_pulse_times;
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

// spin-spin coupling Hamiltonian for NV center with cluster; assumes large static Bz
MatrixXcd H_int_large_static_Bz(const nv_system& nv, const uint cluster_index){
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction between NV center and spin s
    H += act( tp(dot(nv.e.S,zhat), dot(A(nv,cluster.at(s)),nv.nuclei.at(cluster.at(s)).S)),
              {0,s+1}, spins);
    for(uint r = 0; r < s; r++){
      // interaction betwen spin r and spin s
      H += act(H_ss_large_static_Bz(nv.nuclei.at(cluster.at(r)), nv.nuclei.at(cluster.at(s))),
               {r+1,s+1}, spins);
    }
  }
  return H;
}

// NV zero-field splitting plus Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const nv_system& nv, const Vector3d& B){
  return NV_ZFS*dot(nv.e.S,zhat)*dot(nv.e.S,zhat) - nv.e.g*dot(B,nv.e.S);
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
double coherence_measurement(const nv_system& nv, const double w_scan, const uint k_DD,
                             const double f_DD, const double scan_time){
  const double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulse_times = axy_pulse_times(k_DD,f_DD); // AXY protocol pulse times

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

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const uint k_DD, const double f_DD,
                              const double simulation_time, const double axy_pulse_delay){
  // AXY sequence period and pulse_times
  const double t_DD = 2*pi/w_DD;
  const vector<double> pulse_times = delayed_pulse_times(axy_pulse_times(k_DD, f_DD),
                                                         axy_pulse_delay);

  // NV+cluster Hamiltonian
  // const MatrixXcd H = H_int(nv,cluster) + H_Z(nv,cluster,nv.static_Bz*zhat);
  const MatrixXcd H = H_int_large_static_Bz(nv,cluster) + H_nZ(nv,cluster,nv.static_Bz*zhat);

  // NV center spin flip pulse
  const MatrixXcd X = act(sx,{0},nv.clusters.at(cluster).size()+1);

  MatrixXcd U = MatrixXcd::Identity(H.rows(),H.cols());

  // propagator for whole AXY sequences
  if(simulation_time > t_DD){
    for(uint i = 1; i < pulse_times.size(); i++){
      U = (X * exp(-j*H*t_DD*(pulse_times.at(i)-pulse_times.at(i-1))) * U).eval();
    }
    U = pow(X*U, int(simulation_time/t_DD));
  }

  // propagator for AXY sequence remainder
  const double remaining_time = simulation_time - int(simulation_time/t_DD)*t_DD;
  uint pulse_count = 0;
  for(uint i = 1; i < pulse_times.size(); i++){
    if(pulse_times.at(i)*t_DD < remaining_time){
      U = (X * exp(-j*H*t_DD*(pulse_times.at(i)-pulse_times.at(i-1))) * U).eval();
      pulse_count++;
    } else{
      U = (exp(-j*H*(remaining_time-pulse_times.at(i-1)*t_DD)) * U).eval();
      break;
    }
  }
  if(pulse_count%2 != 0) U = (X*U).eval();

  // normalize propagator
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));

  return U;
}

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const uint k_DD, const double f_DD,
                              const double simulation_time, const control_fields& controls,
                              const double axy_pulse_delay){
  const uint spins = nv.clusters.at(cluster).size()+1;

  // AXY sequence period and pulse_times
  const double t_DD = 2*pi/w_DD;
  const vector<double> pulse_times = delayed_pulse_times(axy_pulse_times(k_DD, f_DD),
                                                         axy_pulse_delay);

  // maximim frequency scale of simulation
  double max_freq_scale = w_DD;
  for(uint c = 0; c < controls.num(); c++){
    max_freq_scale = max(max_freq_scale,controls.freqs.at(c));
  }

  // integration step size and number
  const double dt = 1/(max_freq_scale*nv.scale_factor);
  const double dx = 1/(t_DD*max_freq_scale*nv.scale_factor);
  const uint integration_steps = int(simulation_time/dt+0.5);

  // static NV+cluster Hamiltonian
  const MatrixXcd H_static = H_int_large_static_Bz(nv,cluster);

  // initial propagator
  MatrixXcd U = MatrixXcd::Identity(pow(2,spins),pow(2,spins));

  // NV center spin flip pulse
  const MatrixXcd X = act(sx, {0}, spins);

  const uint print_steps = int(nv.scale_factor);
  cout << "Progress (out of " << print_steps << ")...";

  uint pulse_count = 0;
  for(uint x_i = 1; x_i <= integration_steps; x_i++){
    const double x = x_i*dx; // normalized time

    if(x_i%(integration_steps/print_steps) == 0){
      cout << " " << int(print_steps*double(x_i)/integration_steps) << flush;
    }

    // normzlized time into current AXY half-sequence
    const double x_hAXY = x - floor(x/0.5)*0.5; // time in current AXY half-sequence
    // if we are within dx/2 of an AXY pulse time, flip the projections
    if(min({abs(x_hAXY-pulse_times.at(1)), abs(x_hAXY-pulse_times.at(2)),
            abs(x_hAXY-pulse_times.at(3)), abs(x_hAXY-pulse_times.at(4)),
            abs(x_hAXY-pulse_times.at(5))}) < dx/2){
      U = (X*U).eval();
      pulse_count++;
    }

    // net magnetic field
    Vector3d B = nv.static_Bz*zhat;
    for(uint c = 0; c < controls.num(); c++){
      B += controls.Bs.at(c) * cos(controls.freqs.at(c)*x*t_DD + controls.phases.at(c));
    }

    // current Hamiltonian
    const MatrixXcd H = H_static + H_nZ(nv,cluster,B);

    // update propagator
    U = (exp(-j*dx*t_DD*H)*U).eval();
  }
  if(pulse_count%2 != 0) U = (X*U).eval();

  cout << endl << endl;

  // normalize propagator
  U /= sqrt(real(trace(U.adjoint()*U)/double(U.rows())));

  return U;
}

//--------------------------------------------------------------------------------------------
// Single nuclear targeting
//--------------------------------------------------------------------------------------------

// return "natural" basis of a nucleus
vector<Vector3d> natural_basis(const nv_system& nv, const uint index){
  // larmor frequency of and hyperfine field at target nucleus
  const Vector3d larmor_eff = effective_larmor(nv,index);
  const Vector3d hyperfine_perp = A_perp(nv,index);
  // "natural" basis vectors for target nucleus
  const Vector3d target_zhat = hat(larmor_eff);
  const Vector3d target_xhat = hat(hyperfine_perp);
  const Vector3d target_yhat = target_zhat.cross(target_xhat);
  return {target_xhat, target_yhat, target_zhat};
}

// return axis with given polar and azimuthal angles in the "natural" basis
Vector3d natural_axis(const nv_system& nv, const uint index,
                      const double azimuth, const double polar){
  const vector<Vector3d> basis = natural_basis(nv,index);
  return cos(polar) * basis.at(2) + sin(polar) * ( cos(azimuth) * basis.at(0) +
                                                   sin(azimuth) * basis.at(1) );
}

// return control field for decoupling spin s from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd){
  const spin s = nv.nuclei.at(index);
  const Vector3d w_j = effective_larmor(nv,index);
  const double w_rfd = w_j.norm()/(1-sin(theta_rfd)/(2*sqrt(2)*nv.scale_factor));
  const double V_rfd = w_rfd/(s.g*nv.scale_factor);
  const Vector3d n_rfd = cos(theta_rfd)*hat(w_j) + sin(theta_rfd)*hat(w_j.cross(zhat));
  return control_fields(V_rfd*n_rfd, w_rfd, phi_rfd);
}

// exact gate G = exp(-i * rotation_angle * sigma_{axis}^{index})
MatrixXcd G_ctl(const nv_system& nv, const uint index, const double target_axis_azimuth,
                const double rotation_angle){
  const uint spins = nv.clusters.at(get_cluster_containing_index(nv,index)).size()+1;
  const uint index_in_cluster = get_index_in_cluster(nv,index);
  const Vector3d target_axis = natural_axis(nv, index, target_axis_azimuth);
  const MatrixXcd G = exp(-j * rotation_angle * dot(s_vec,target_axis));
  return act(G, {index_in_cluster+1}, spins);
}

// approximate propagator U = exp(-i * rotation_angle * sigma_{axis}^{index})
MatrixXcd U_ctl(const nv_system& nv, const uint index, const double target_axis_azimuth,
                const double rotation_angle){
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_index(nv, index);
  const uint spins = nv.clusters.at(cluster).size()+1;

  for(uint i = 0; i < nv.clusters.at(cluster).size(); i++){
    if(nv.clusters.at(cluster).at(i) == index) continue;
    if(is_larmor_pair(nv, index, nv.clusters.at(cluster).at(i))){
      cout << "Cannot address nuclei with larmor pairs: "
           << index << ", " << nv.clusters.at(cluster).at(i) << endl;
      return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
    }
  }

  // larmor frequency of target nucleus
  const double w_larmor = effective_larmor(nv,index).norm();
  const double t_larmor = 2*pi/w_larmor;

  // axis of rotation
  const Vector3d axis = natural_axis(nv, index, target_axis_azimuth);

  // AXY protocol parameters
  const double A_j = A(nv,index).norm();
  double w_DD;
  if(w_larmor >= nv.scale_factor*A_j){
    uint k_M = int(w_larmor/(nv.scale_factor*A_j)-1);
    if(k_M % 2 == 0) k_M--;
    w_DD = (w_larmor-nv.scale_factor*A_j)/k_M;
  }
  if(w_larmor < nv.scale_factor*A_j || w_DD < nv.scale_factor*A_j){
    w_DD = (w_larmor+nv.scale_factor*A_j)/3.;
  }
  const double t_DD = 2*pi/w_DD;
  const uint k_DD = abs(w_DD - w_larmor) < abs(3*w_DD - w_larmor) ? 1 : 3;
  const double f_DD = 0;

  // control field frequency = effective larmor frequency of target nucleus
  const double w_ctl = w_larmor; // control field frequency
  const double dw_min = larmor_resolution(nv,index);
  double g_B_ctl = dw_min/nv.scale_factor; // ctl field strength * gyromangnetic ratio

  const double control_period = 2*pi*4/g_B_ctl;
  double control_time = -4*rotation_angle/g_B_ctl; // control operation time
  while(control_time >= control_period) control_time -= control_period;
  while(control_time < 0) control_time += control_period;
  if(control_time > control_period/2){
    g_B_ctl *= -1;
    control_time = control_period-control_time;
  }
  const double flush_time = ceil(control_time/t_larmor)*t_larmor - control_time;

  const double B_ctl = g_B_ctl/nv.nuclei.at(index).g; // control field strength
  const control_fields controls(B_ctl*axis, w_ctl, 0.); // control field object

  const MatrixXcd U_control = simulate_propagator(nv, cluster, w_DD, k_DD, f_DD,
                                                  control_time, controls);

  const double axy_delay = control_time/t_DD - int(control_time/t_DD);
  const MatrixXcd U_flush = simulate_propagator(nv, cluster, w_DD, k_DD, f_DD,
                                                flush_time, axy_delay);

  return U_flush * U_control;
}

// exact gate G = exp(-i * rotation_angle * sigma_{n_1}^{NV}*sigma_{n_2}^{index})
MatrixXcd G_int(const nv_system& nv, const uint index, const uint k_DD,
                const double nv_axis_polar, const double nv_axis_azimuth,
                const double target_axis_azimuth, const double rotation_angle){
  // spin of NV center along given axis
  const Vector3d nv_axis = natural_axis(nv, index, nv_axis_azimuth, nv_axis_polar);
  const Matrix2cd sn_nv = dot(s_vec,nv_axis);

  // spin of target nucleus along given axis
  const vector<Vector3d> target_basis = natural_basis(nv,index);
  const Matrix2cd sxp = dot(s_vec,target_basis.at(0));
  const Matrix2cd syp = dot(s_vec,target_basis.at(1));
  const Matrix2cd sn_target = cos(target_axis_azimuth)*sxp + sin(target_axis_azimuth)*syp;

  // exact coupling operation
  const uint spins = nv.clusters.at(get_cluster_containing_index(nv,index)).size()+1;
  const uint index_in_cluster = get_index_in_cluster(nv,index);
  const MatrixXcd G = exp(-j * rotation_angle * tp(sn_nv,sn_target));
  return act(G, {0,index_in_cluster+1}, spins);
}

// approximate propagator U = exp(-i * rotation_angle * sigma_{n_1}^{NV}*sigma_{n_2}^{index})
MatrixXcd U_int(const nv_system& nv, const uint index, const uint k_DD,
                const double nv_axis_polar, const double nv_axis_azimuth,
                const double target_axis_azimuth, const double rotation_angle){
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_index(nv, index);
  const uint index_in_cluster = get_index_in_cluster(nv, index);
  const uint spins = nv.clusters.at(cluster).size()+1;

  // verify that we can address this nucleus
  if(round(4*dot(nv.nuclei.at(index).pos,ao)) == 0.){
    cout << "Cannot address nuclei without hyperfine coupling perpendicular to the NV axis: "
         << index << endl;
    return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  }
  for(uint i = 0; i < nv.clusters.at(cluster).size(); i++){
    if(nv.clusters.at(cluster).at(i) == index) continue;
    if(is_larmor_pair(nv, index, nv.clusters.at(cluster).at(i))){
      cout << "Cannot address nuclei with larmor pairs: "
           << index << ", " << nv.clusters.at(cluster).at(i) << endl;
      return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
    }
  }

  // larmor frequency of and perpendicular component of hyperfine field at target nucleus
  const double w_larmor = effective_larmor(nv,index).norm();
  const double dw_min = larmor_resolution(nv,index);
  const Vector3d hyperfine_perp = A_perp(nv,index);

  // AXY sequence parameters
  const double w_DD = w_larmor/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  double f_DD = dw_min/(hyperfine_perp.norm()*nv.scale_factor);

  const double interaction_period = 2*pi/abs(f_DD*hyperfine_perp.norm()/8);
  double interaction_time = rotation_angle/(nv.ms*f_DD*hyperfine_perp.norm()/8);
  while(interaction_time >= interaction_period) interaction_time -= interaction_period;
  while(interaction_time < 0) interaction_time += interaction_period;
  if(interaction_time > interaction_period/2){
    f_DD *= -1;
    interaction_time = interaction_period-interaction_time;
  }
  const double flush_time = ceil(interaction_time/t_DD)*t_DD - interaction_time;

  const MatrixXcd U_interaction =
    simulate_propagator(nv, cluster, w_DD, k_DD, f_DD,
                        interaction_time, -target_axis_azimuth/(2*pi));

  const double axy_delay = interaction_time/t_DD - int(interaction_time/t_DD);
  const MatrixXcd U_flush =
    simulate_propagator(nv, cluster, w_DD, k_DD, 0,
                        flush_time, axy_delay - target_axis_azimuth/(2*pi));

  const Vector3d nv_axis = natural_axis(nv, index, nv_axis_azimuth, nv_axis_polar);

  // rotate the NV spin between the desired axis and zhat
  const MatrixXcd nv_axis_to_zhat = act( rotate(zhat, nv_axis), {0}, spins);

  return
    U_flush *
    nv_axis_to_zhat.adjoint() *
    U_interaction *
    nv_axis_to_zhat;
}

// compute fidelity of iSWAP operation between NV center and target nucleus
double iSWAP_fidelity(const nv_system& nv, const uint index, const uint k_DD){
  const double iswap_angle = -pi/4;
  const double nv_axis_polar = pi/2;
  const double xhat_azimuth = 0;
  const double yhat_azimuth = pi/2;

  // approximate iSWAP gate
  const MatrixXcd U_sx = U_int(nv, index, k_DD, nv_axis_polar,
                               xhat_azimuth, xhat_azimuth, iswap_angle);
  const MatrixXcd U_sy = U_int(nv, index, k_DD, nv_axis_polar,
                               yhat_azimuth, yhat_azimuth, iswap_angle);
  const MatrixXcd U_iSWAP = U_sy * U_sx;

  // exact iSWAP gape
  const MatrixXcd G_sx = G_int(nv, index, k_DD, nv_axis_polar,
                               xhat_azimuth, xhat_azimuth, iswap_angle);
  const MatrixXcd G_sy = G_int(nv, index, k_DD, nv_axis_polar,
                               yhat_azimuth, yhat_azimuth, iswap_angle);
  const MatrixXcd iSWAP = G_sy * G_sx;

  return gate_fidelity(U_iSWAP, iSWAP);
}

// return SWAP operation between NV center and the ST subspace of two nuclei
double SWAP_NVST_fidelity(const nv_system& nv, const uint idx1, const uint idx2,
                          const uint k_DD){
  // assert that both nuclei are in the same
  const uint cluster = get_cluster_containing_index(nv,idx1);
  assert(in_vector(idx2,nv.clusters.at(cluster)));

  const uint spins = nv.clusters.at(cluster).size()+1;

  const Matrix2cd R_n1 = rotate(natural_basis(nv,idx1),{xhat,yhat,zhat});
  const Matrix2cd R_n2 = rotate(natural_basis(nv,idx2),{xhat,yhat,zhat});
  const MatrixXcd NV_to_n1 = act(R_n1,{0},spins);
  const MatrixXcd NV_to_n2 = act(R_n2,{0},spins);

  const MatrixXcd Rz_NV = act(U_NV(zhat,pi/4),{0},spins);

  // exact SWAP_NVST gate
  const MatrixXcd Rx_1_exact = NV_to_n1.adjoint() * G_ctl(nv,idx1,0,pi/4) * NV_to_n1;
  const MatrixXcd Ry_1_exact = NV_to_n1.adjoint() * G_ctl(nv,idx1,pi/2,pi/4) * NV_to_n1;
  const MatrixXcd Rz_1_exact = Rx_1_exact * Ry_1_exact * Rx_1_exact.adjoint();
  const MatrixXcd iSWAP_NV_1_exact =
    NV_to_n1.adjoint() * G_int(nv,idx1,k_DD,pi/2,0,0,-pi/4) *
    G_int(nv,idx1,k_DD,pi/2,pi/2,pi/2,-pi/4) * NV_to_n1;
  const MatrixXcd E_NV_2_exact =
    NV_to_n2.adjoint() * G_int(nv,idx2,k_DD,pi/2,pi/2,0,-pi/4) * NV_to_n2;
  const MatrixXcd cNOT_NV_1_exact =
    Rz_NV * NV_to_n1.adjoint() * Rx_1_exact * G_int(nv,idx1,k_DD,0,0,0,-pi/4) * NV_to_n1;

  const MatrixXcd SWAP_NVST =
    Rz_NV * Rz_1_exact * iSWAP_NV_1_exact.adjoint() *
    E_NV_2_exact * cNOT_NV_1_exact * E_NV_2_exact.adjoint() *
    iSWAP_NV_1_exact * Rz_1_exact.adjoint() * Rz_NV.adjoint();

  // approximate SWAP_NVST gate
  const MatrixXcd Rx_1 = NV_to_n1.adjoint() * U_ctl(nv,idx1,0,pi/4) * NV_to_n1;
  const MatrixXcd Ry_1 = NV_to_n1.adjoint() * U_ctl(nv,idx1,pi/2,pi/4) * NV_to_n1;
  const MatrixXcd Rz_1 = Rx_1 * Ry_1 * Rx_1.adjoint();
  const MatrixXcd iSWAP_NV_1 =
    NV_to_n1.adjoint() * U_int(nv,idx1,k_DD,pi/2,0,0,-pi/4) *
    U_int(nv,idx1,k_DD,pi/2,pi/2,pi/2,-pi/4) * NV_to_n1;
  const MatrixXcd E_NV_2 =
    NV_to_n2.adjoint() * U_int(nv,idx2,k_DD,pi/2,pi/2,0,-pi/4) * NV_to_n2;
  const MatrixXcd cNOT_NV_1 =
    Rz_NV * NV_to_n1.adjoint() * Rx_1 * U_int(nv,idx1,k_DD,0,0,0,-pi/4) * NV_to_n1;

  const MatrixXcd U_SWAP_NVST =
    Rz_NV * Rz_1 * iSWAP_NV_1.adjoint() *
    E_NV_2 * cNOT_NV_1 * E_NV_2.adjoint() *
    iSWAP_NV_1 * Rz_1.adjoint() * Rz_NV.adjoint();

  return gate_fidelity(U_SWAP_NVST,SWAP_NVST);
}
