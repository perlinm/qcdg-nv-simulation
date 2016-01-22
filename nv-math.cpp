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

nv_system::nv_system(const int ms, const double static_Bz, const double scale_factor) :
  e(spin(Vector3d::Zero(), ge,
         mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat))),
  ms(ms), static_Bz(static_Bz), scale_factor(scale_factor)
{};

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// determine whether two spins are a larmor pair
bool larmor_pair(const spin& s1, const spin& s2, const double tolerance){

  const double r1z = dot(s1.pos,zhat);
  const double r2z = dot(s2.pos,zhat);

  const double r1xy = (s1.pos - r1z*zhat).norm();
  const double r2xy = (s2.pos - r2z*zhat).norm();

  return abs(abs(r1z/r2z)-1) < tolerance && abs(abs(r1xy/r2xy)-1) < tolerance;
}

// determine whether two spins are in the same larmor group
bool larmor_group(const nv_system& nv, const uint idx1, const uint idx2){
  const double dw = (effective_larmor(nv,idx2).norm() - effective_larmor(nv,idx1).norm());
  return abs(dw) < nv.cluster_coupling/nv.scale_factor;
}

// coupling strength between two spins; assumes strong magnetic field in zhat
inline double coupling_strength(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos - s1.pos;
  return s1.g*s2.g/(8*pi*pow(r.norm()*a0,3)) * (1 - 3*dot(hat(r),zhat)*dot(hat(r),zhat));
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
          if(larmor_group(nv, new_cluster.at(i), old_cluster.at(j))){
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

// compute fidelity of iSWAP operation between NV center and target nucleus
fidelity_info iswap_fidelity(const nv_system& nv, const uint index, const uint k_DD){
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_index(nv, index);
  const uint index_in_cluster = get_index_in_cluster(nv, index);
  const uint spins = nv.clusters.at(cluster).size()+1;

  // larmor frequency of and hyperfine field at target nucleus
  const Vector3d larmor_eff = effective_larmor(nv,index);
  const Vector3d hyperfine = A(nv,index);
  const Vector3d hyperfine_w = dot(hyperfine,hat(larmor_eff))*hat(larmor_eff);
  const Vector3d hyperfine_perp = hyperfine - hyperfine_w;

  // if our interaction strength is too weak, we can't address this nucleus
  if(hyperfine_perp.norm() < nv.cluster_coupling) return fidelity_info();

  // minimum difference in larmor frequencies between target nucleus and other nuclei
  double dw_min = DBL_MAX;
  for(uint s = 0; s < nv.nuclei.size(); s++){
    if(s == index) continue;
    const double dw = abs(larmor_eff.norm() - effective_larmor(nv,s).norm());
    if(dw < dw_min) dw_min = dw;
  }

  // if the target's larmor frequency is too close to another,
  //   it is unaddressable via the default protocol
  if(dw_min < nv.cluster_coupling/nv.scale_factor) return fidelity_info();

  // AXY sequence parameters
  const double f_DD = -nv.ms*dw_min/(hyperfine_perp.norm()*nv.scale_factor);
  const double w_DD = larmor_eff.norm()/k_DD; // AXY protocol angular frequency
  const double operation_time = 2*pi/abs(f_DD*hyperfine_perp.norm());

  cout << "Running iSWAP protocol with the following parameters:"
       << " target nucleus index: " << index << endl
       << " w_DD: " << w_DD << endl
       << " f_" << k_DD << ": " << f_DD << endl
       << " operation_time: " << operation_time*1000 << " ms" << endl;

  // spin addressing propagators

  const MatrixXcd U_sx = simulate_propagator(nv, cluster, w_DD, k_DD, f_DD, operation_time);
  const MatrixXcd U_sy = simulate_propagator(nv, cluster, w_DD, k_DD, f_DD,
                                             operation_time, 0.25);

  // "natural" basis vectors for target nucleus
  const Vector3d target_zhat = hat(larmor_eff);
  const Vector3d target_xhat = hat(hyperfine_perp);
  const Vector3d target_yhat = target_zhat.cross(target_xhat);

  // operators to rotate into the "natural" frame
  const Matrix2cd Z_to_tY = rotate(acos(dot(zhat,target_yhat)), zhat.cross(target_yhat));
  const Matrix2cd Z_to_tX = rotate(acos(dot(zhat,target_xhat)), zhat.cross(target_xhat));

  // propagator for the iSWAP operation
  const MatrixXcd U_iSWAP =
    act(Z_to_tX,{0},spins) * U_sx * act(Z_to_tX.inverse(),{0},spins) *
    act(Z_to_tY,{0},spins) * U_sy * act(Z_to_tY.inverse(),{0},spins);

  // spin operators in the "natural" basis
  const Matrix2cd sxp = dot(s_vec,target_xhat);
  const Matrix2cd syp = dot(s_vec,target_yhat);

  // exact iSWAP gate
  const MatrixXcd iSWAP = exp(j*0.25*pi * act( tp(sxp,sxp) + tp(syp,syp),
                                               {0,index_in_cluster+1}, spins));

  // compute iSWAP gate fidelity
  const double fidelity = gate_fidelity(U_iSWAP,iSWAP);

  return fidelity_info(larmor_eff.norm(), hyperfine.norm(), hyperfine_perp.norm(),
                       dw_min, f_DD, operation_time, fidelity);
}

// realization of propagator U = exp(-i * phi * sigma_{axis}^{index})
MatrixXcd U_ctl(const nv_system& nv, const uint index, const double target_axis_azimuth,
                const double phi, const double threshold){
  const uint cluster = get_cluster_containing_index(nv, index);
  const uint spins = nv.clusters.at(cluster).size()+1;

  // larmor frequency of and hyperfine field at target nucleus
  const Vector3d larmor_eff = effective_larmor(nv,index);
  const Vector3d hyperfine = A(nv,index);
  const Vector3d hyperfine_perp = hyperfine - dot(hyperfine,hat(larmor_eff))*hat(larmor_eff);

  const double w_larmor = larmor_eff.norm();
  const double t_larmor = 2*pi/w_larmor;

  // "natural" basis vectors for target nucleus
  const Vector3d target_zhat = hat(larmor_eff);
  const Vector3d target_xhat = hat(hyperfine_perp);
  const Vector3d target_yhat = target_zhat.cross(target_xhat);

  const Vector3d axis =
    cos(target_axis_azimuth)*target_xhat + sin(target_axis_azimuth)*target_yhat;

  // control field frequency = effective larmor frequency of target nucleus

  double dw_min = DBL_MAX; // min{ |w_s - w_{target}| for all s != index }
  for(uint s = 0; s < nv.nuclei.size(); s++){
    if(s == index) continue;
    const double dw = abs(w_larmor - effective_larmor(nv,s).norm());
    if(dw < dw_min) dw_min = dw;
  }

  if(dw_min < threshold){
    cout << "Cannot target nuclei with larmor pairs"
         << " (target nucleus: " << index << ")" << endl;
    return MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  }

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

  const double w_ctl = w_larmor; // control field frequency
  double g_B_ctl = dw_min/nv.scale_factor; // ctl field strength * gyromangnetic ratio

  const double control_period = 2*pi*4/g_B_ctl;
  double control_time = 4*phi/g_B_ctl; // control operation time
  while(control_time >= control_period) control_time -= control_period;
  while(control_time < 0) control_time += control_period;
  if(control_time > control_period/2){
    g_B_ctl *= -1;
    control_time = control_period-control_time;
  }
  const double flush_time = ceil(control_time/t_larmor)*t_larmor - control_time;

  const double B_ctl = g_B_ctl/nv.nuclei.at(index).g; // control field strength
  const control_fields controls(B_ctl*axis, w_ctl, 0.); // control field object

  cout << "Running targeting protocol with the following parameters:" << endl
       << " target nucleus index: " << index << endl
       << " w_DD: " << w_DD << endl
       << " f_" << k_DD << ": 0" << endl
       << " w_ctl: " << w_ctl << endl
       << " B_ctl: " << B_ctl/gauss*1000 << " mG" << endl
       << " control_time: " << control_time*1000 << " ms" << endl
       << endl;

  cout << "Additional information:" << endl
       << " dw_min: " << dw_min << endl
       << " A: " << A_j << endl
       << " cluster size: " << nv.clusters.at(cluster).size() << endl
       << endl;

  const MatrixXcd U_control = simulate_propagator(nv, cluster, w_DD, k_DD, f_DD,
                                                  control_time, controls);

  const double axy_delay = control_time/t_DD - int(control_time/t_DD);
  const MatrixXcd U_flush = simulate_propagator(nv, cluster, w_DD, k_DD, f_DD,
                                                flush_time, axy_delay);

  return U_flush * U_control;
}

// propagator generated by t_0 H = phi * sigma_{nv_axis}^{NV} * sigma_{axis_ind}^{index}
MatrixXcd U_int(const nv_system& nv, const uint index, const uint k_DD,
                const Vector3d& nv_axis, const double target_axis_azimuth, const double phi,
                const double threshold){

  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_index(nv, index);
  const uint index_in_cluster = get_index_in_cluster(nv, index);
  const uint spins = nv.clusters.at(cluster).size()+1;

  // larmor frequency of and hyperfine field at target nucleus
  const Vector3d larmor_eff = effective_larmor(nv,index);
  const Vector3d hyperfine = A(nv,index);
  const Vector3d hyperfine_perp = hyperfine - dot(hyperfine,hat(larmor_eff))*hat(larmor_eff);

  const double w_larmor = larmor_eff.norm();

  // "natural" basis vectors for target nucleus
  const Vector3d target_zhat = hat(larmor_eff);
  const Vector3d target_xhat = hat(hyperfine_perp);
  const Vector3d target_yhat = target_zhat.cross(target_xhat);

  // minimum difference in larmor frequencies between target nucleus and other nuclei
  double dw_min = DBL_MAX;
  for(uint s = 0; s < nv.nuclei.size(); s++){
    if(s == index) continue;
    const double dw = abs(w_larmor - effective_larmor(nv,s).norm());
    if(dw < dw_min) dw_min = dw;
  }

  // AXY sequence parameters
  const double w_DD = w_larmor/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  double f_DD = dw_min/(hyperfine_perp.norm()*nv.scale_factor);

  const double interaction_period = 2*pi/abs(f_DD*hyperfine_perp.norm()/8);
  double interaction_time = phi/(nv.ms*f_DD*hyperfine_perp.norm()/8);
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

  const Vector3d target_nv_axis =
    dot(nv_axis,xhat)*target_xhat +
    dot(nv_axis,yhat)*target_yhat +
    dot(nv_axis,zhat)*target_zhat;

  // operator to rotate into the "natural" frame
  const Matrix2cd rotate_to_target_frame =
    rotate(acos(dot(zhat,target_nv_axis)),zhat.cross(target_nv_axis));

  // return
    // act(rotate_to_target_frame.inverse(),{0},spins) *
    // U_flush * U_interaction *
    // act(rotate_to_target_frame,{0},spins);

  const MatrixXcd U =
    act(rotate_to_target_frame.inverse(),{0},spins) *
    U_flush * U_interaction *
    act(rotate_to_target_frame,{0},spins);

  // spin operators in the "natural" basis
  const Matrix2cd sxp = dot(s_vec,target_xhat);
  const Matrix2cd syp = dot(s_vec,target_yhat);
  const Matrix2cd szp = dot(s_vec,target_zhat);

  const Matrix2cd sn_target = cos(target_axis_azimuth)*sxp + sin(target_axis_azimuth)*syp;
  const Matrix2cd sn_nz =
    dot(nv_axis,xhat)*sxp +
    dot(nv_axis,yhat)*syp +
    dot(nv_axis,zhat)*szp;

  // desired propagator
  const MatrixXcd U_desired =
    exp(-j*phi * act(tp(sn_nz,sn_target), {0,index_in_cluster+1}, spins));

  cout << remove_artifacts(remove_phase(U),1e-3) << endl << endl;
  cout << remove_artifacts(remove_phase(U_desired),1e-3) << endl << endl;
  cout << gate_fidelity(U,U_desired) << endl << endl;

  return U;
}

