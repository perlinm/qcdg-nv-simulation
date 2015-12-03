#include <iostream> // for standard output
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "constants.h"
#include "qp-math.h"
#include "nv-math.h"

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// coupling strength between two spins; assumes strong magnetic field in zhat
inline double coupling_strength(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos - s1.pos;
  return s1.g*s2.g/(8*pi*pow(r.norm()*a0,3)) * (1 - 3*dot(hat(r),zhat)*dot(hat(r),zhat));
}

// group nuclei into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<uint>> get_index_clusters(const vector<spin>& nuclei,
                                        const double min_coupling_strength){

  vector<vector<uint>> ind_clusters; // all clusters
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
      ind_clusters.push_back(cluster);
    }
  }
  return ind_clusters;
}

// get size of largest spin cluster
uint largest_cluster_size(const vector<vector<uint>>& ind_clusters){
  uint largest_size = 0;
  for(uint i = 0; i < ind_clusters.size(); i++){
    if(ind_clusters.at(i).size() > largest_size) largest_size = ind_clusters.at(i).size();
  }
  return largest_size;
}

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(const vector<spin>& nuclei, const uint cluster_size_target,
                            const double initial_cluster_coupling, const double dcc_cutoff){

  double cluster_coupling = initial_cluster_coupling;
  double dcc = cluster_coupling/4;

  vector<vector<uint>> ind_clusters = get_index_clusters(nuclei, cluster_coupling);
  bool coupling_too_small = largest_cluster_size(ind_clusters) >= cluster_size_target;
  bool last_coupling_too_small;
  bool crossed_correct_coupling;

  // determine coupling for which the largest cluster size just barely >= max_cluster_size
  while(dcc >= dcc_cutoff || !coupling_too_small){
    last_coupling_too_small = coupling_too_small;

    cluster_coupling += coupling_too_small ? dcc : -dcc;
    ind_clusters = get_index_clusters(nuclei,cluster_coupling);
    coupling_too_small = largest_cluster_size(ind_clusters) >= cluster_size_target;

    if(coupling_too_small != last_coupling_too_small) crossed_correct_coupling = true;

    if(coupling_too_small == last_coupling_too_small){
      if(!crossed_correct_coupling) dcc *= 2;
    } else{
      dcc /= 2;
    }
  }
  return cluster_coupling;
}

// group nuclei into clusters according to cluster_indices
vector<vector<spin>> group_nuclei(const vector<spin>& nuclei,
                                 const vector<vector<uint>>& cluster_indices){
  vector<vector<spin>> clusters;
  for(uint c = 0; c < cluster_indices.size(); c++){
    vector<spin> cluster;
    for(uint s = 0; s < cluster_indices.at(c).size(); s++){
      cluster.push_back(nuclei.at(cluster_indices.at(c).at(s)));
    }
    clusters.push_back(cluster);
  }
  return clusters;
}

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// pulse times for harmonic k, fourier component f, and a given offset
vector<double> axy_pulses(const uint k, const double f, double offset = 0.){

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
  vector<double> pulses;
  pulses.push_back(x1);
  pulses.push_back(x2);
  pulses.push_back(0.25);
  pulses.push_back(0.5 - x2);
  pulses.push_back(0.5 - x1);
  pulses.push_back(0.5 + x1);
  pulses.push_back(0.5 + x2);
  pulses.push_back(0.75);
  pulses.push_back(1. - x2);
  pulses.push_back(1. - x1);

  if(offset == 0.) return pulses;

  // shift offset to be within one AXY period
  while(offset < 0) offset += 1;
  while(offset > 1) offset -= 1;

  vector<double> offset_pulses;
  for(uint p = 0; p < 20; p++){
    if(p/10 + pulses.at(p%10) - offset > 0){
      offset_pulses.push_back(p/10 + pulses.at(p%10) - offset);
    }
    if(offset_pulses.size() == 10) break;
  }
  return offset_pulses;
}

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2){
  const Vector3d r = s2.pos-s1.pos;
  return s1.g*s2.g/(4*pi*pow(r.norm()*a0,3))
    * (dot(s1.S,s2.S) - 3*tp(dot(s1.S,hat(r)), dot(s2.S,hat(r))));
}

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const spin& e, const vector<spin>& cluster){
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction between NV center and spin s
    H += act(H_ss(e, cluster.at(s)), {0,s+1}, spins);
    for(uint r = 0; r < s; r++){
      // interaction betwen spin r and spin s
      H += act(H_ss(cluster.at(r), cluster.at(s)), {r+1,s+1}, spins);
    }
  }
  return H;
}

// NV zero-field splitting plus Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const spin& e, const vector<spin>& cluster, const Vector3d& B){
  return act(NV_ZFS*dot(e.S,zhat)*dot(e.S,zhat) - e.g*dot(B,e.S), {0}, cluster.size()+1);
}

// nuclear Zeeman Hamiltonian
MatrixXcd H_nZ(const vector<spin>& cluster, const Vector3d& B){
  const int spins = cluster.size()+1;
  // zero-field splitting and interaction of NV center with magnetic field
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction of spin s with magnetic field
    H -= act(cluster.at(s).g*dot(B,cluster.at(s).S), {s+1}, spins);
  }
  return H;
}

// Zeeman Hamiltonian for NV+cluster system
inline MatrixXcd H_Z(const spin& e, const vector<spin>& cluster, const Vector3d& B){
  return H_nZ(cluster, B) + H_NV_GS(e, cluster, B);
}

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const int ms, const vector<vector<spin>>& clusters,
                             const double w_scan, const uint k_DD, const double f_DD,
                             const double scan_time, const double static_B){
  const spin e_ms = e(ms); // electron spin

  const double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulses = axy_pulses(k_DD,f_DD); // AXY protocol pulse times

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    // projections onto |ms> and |0> NV states
    const MatrixXcd proj_m = act(up*up.adjoint(), {0}, clusters.at(c).size()+1); // |ms><ms|
    const MatrixXcd proj_0 = act(dn*dn.adjoint(), {0}, clusters.at(c).size()+1); // |0><0|

    // construct full Hamiltonian
    const MatrixXcd H = H_int(e_ms,clusters.at(c)) + H_Z(e_ms,clusters.at(c),static_B*zhat);

      // spin cluster Hamiltonians for single NV states
    const MatrixXcd H_m = ptrace(H*proj_m, {0}); // <ms|H|ms>
    const MatrixXcd H_0 = ptrace(H*proj_0, {0}); // <0|H|0>

    // propagators for sections of the AXY sequence
    const MatrixXcd U1_m = exp(-j*H_m*t_DD*pulses.at(0));
    const MatrixXcd U2_m = exp(-j*H_0*t_DD*(pulses.at(1)-pulses.at(0)));
    const MatrixXcd U3_m = exp(-j*H_m*t_DD*(pulses.at(2)-pulses.at(1)));

    const MatrixXcd U1_0 = exp(-j*H_0*t_DD*pulses.at(0));
    const MatrixXcd U2_0 = exp(-j*H_m*t_DD*(pulses.at(1)-pulses.at(0)));
    const MatrixXcd U3_0 = exp(-j*H_0*t_DD*(pulses.at(2)-pulses.at(1)));

    // single AXY sequence propagators
    MatrixXcd U_m = U1_m*U2_m*U3_m*U3_0*U2_0*U1_0 * U1_0*U2_0*U3_0*U3_m*U2_m*U1_m;
    MatrixXcd U_0 = U1_0*U2_0*U3_0*U3_m*U2_m*U1_m * U1_m*U2_m*U3_m*U3_0*U2_0*U1_0;

    // propagators for entire scan
    U_m = pow(U_m, int(scan_time/t_DD));
    U_0 = pow(U_0, int(scan_time/t_DD));

    // update coherence
    coherence *= real(trace(U_0.adjoint()*U_m)) / pow(2,clusters.at(c).size());
  }
  return coherence;
}

//--------------------------------------------------------------------------------------------
// Control field scanning
// (assume large static magnetic field along NV axis)
//--------------------------------------------------------------------------------------------

// Hamiltoninan coupling two spins
inline MatrixXcd H_ss_large_static_B(const spin& s1, const spin& s2){
  return coupling_strength(s1,s2) * (3*tp(dot(s1.S,zhat), dot(s2.S,zhat)) - dot(s1.S,s2.S));
}

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int_large_static_B(const spin& e, const vector<spin>& cluster){
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction between NV center and spin s
    H += act( tp(dot(e.S,zhat), dot(A(cluster.at(s)),cluster.at(s).S)), {0,s+1}, spins);
    for(uint r = 0; r < s; r++){
      // interaction betwen spin r and spin s
      H += act(H_ss_large_static_B(cluster.at(r), cluster.at(s)), {r+1,s+1}, spins);
    }
  }
  return H;
}

// perform NV coherence measurement with a static magnetic field and additional control fields
double coherence_measurement(const int ms, const vector<vector<spin>>& clusters,
                             const double w_scan, const uint k_DD, const double f_DD,
                             double scan_time, const double static_B,
                             const control_fields& controls, const uint integration_factor){
  const spin e_ms = e(ms); // electron spin

  const double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulses = axy_pulses(k_DD,f_DD); // AXY protocol pulse times

  Vector3d max_possible_B = static_B*zhat;
  for(uint c = 0; c < controls.num(); c++){
    Vector3d B = controls.Bs.at(c);
    max_possible_B += abs(dot(B,xhat))*xhat + abs(dot(B,yhat))*yhat + abs(dot(B,zhat))*zhat;
  }

  // todo: compute freq_scale, dt, dx, and integration_steps for each cluster
  //   for the event that we simulate heteronuclear systems
  const double freq_scale = max(w_DD,gC13*max_possible_B.norm());
  const double dt = 1/(freq_scale*integration_factor); // integration time step
  const double dx = 1/(t_DD*freq_scale*integration_factor);

  scan_time = int(scan_time/t_DD)*t_DD; // scan an integer number of AXY sequences
  const uint integration_steps = int(scan_time/dt);

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    const int spins = clusters.at(c).size()+1; // number of spins in NV+cluster
    const double cHD = pow(2,spins-1); // dimensionality of cluster Hilbert space

    // projections onto |ms> and |0> NV states
    MatrixXcd proj_m = act(up*up.adjoint(), {0}, spins); // |ms><ms|
    MatrixXcd proj_0 = act(dn*dn.adjoint(), {0}, spins); // |0><0|

    // initialize spin cluster propagators
    const MatrixXcd Id = MatrixXcd::Identity(cHD,cHD); // identity matrix
    MatrixXcd U_m = Id; // <ms|U|ms> (initial)
    MatrixXcd U_0 = Id; // <0|U|0> (initial)

    for(uint x_i = 1; x_i <= integration_steps; x_i++){
      const double x = x_i*dx; // time

      // flip electron spin according to the AXY sequence
      const double x_hAXY = x - floor(x/0.5)*0.5; // time in current AXY half-sequence
      // if we are within dx/2 of an AXY pulse time, flip the projections
      if(min({abs(x_hAXY-pulses.at(0)), abs(x_hAXY-pulses.at(1)), abs(x_hAXY-pulses.at(2)),
              abs(x_hAXY-pulses.at(3)), abs(x_hAXY-pulses.at(4))}) < dx/2){
        const MatrixXcd proj_tmp = proj_m;
        proj_m = proj_0;
        proj_0 = proj_tmp;
      }

      // compute net magnetic field
      Vector3d B = static_B*zhat;
      for(uint c = 0; c < controls.num(); c++){
        B += controls.Bs.at(c) * cos(controls.freqs.at(c)*x*t_DD + controls.phases.at(c));
      }

      // compute NV+cluster Hamiltonian
      const MatrixXcd H = H_int_large_static_B(e_ms, clusters.at(c))
        + H_nZ(clusters.at(c), B);

      // spin cluster Hamiltonians
      const MatrixXcd H_m = ptrace(H*proj_m, {0}); // <ms|H|ms>
      const MatrixXcd H_0 = ptrace(H*proj_0, {0}); // <0|H|0>

      // update and normalize propagators
      // U_m = (exp(-j*dt*H_m)*U_m).eval();
      // U_0 = (exp(-j*dt*H_0)*U_0).eval();
      // U_m = ((Id - j*dt*H_m - dt*dt*H_m*H_m/2)*U_m).eval();
      // U_0 = ((Id - j*dt*H_0 - dt*dt*H_0*H_0/2)*U_0).eval();
      U_m = ((Id - j*dt*H_m)*U_m).eval();
      U_0 = ((Id - j*dt*H_0)*U_0).eval();

      U_m /= sqrt(real(trace(U_m.adjoint()*U_m)/cHD));
      U_0 /= sqrt(real(trace(U_0.adjoint()*U_0)/cHD));
    }

    // update coherence
    coherence *= real(trace(U_0.adjoint()*U_m)) / cHD;
  }
  return coherence;
}

//--------------------------------------------------------------------------------------------
// Single nuclear targeting
//--------------------------------------------------------------------------------------------

// return control field for decoupling spin s from other nuclei
control_fields nuclear_decoupling_field(const spin& s, const double static_B, const int ms,
                                        const double phi_rfd, const double theta_rfd,
                                        const double scale){
  const Vector3d w_j = effective_larmor(s, static_B, ms);
  const double w_rfd = w_j.norm()/(1-sin(theta_rfd)/(2*sqrt(2)*scale));
  const double V_rfd = w_rfd/(s.g*scale);
  const Vector3d n_rfd = cos(theta_rfd)*hat(w_j) + sin(theta_rfd)*hat(w_j.cross(zhat));
  return control_fields(V_rfd*n_rfd, w_rfd, phi_rfd);
}

// compute fidelity of SWAP operation between NV center and target spin
double iswap_fidelity(const MatrixXcd& rho_NV_0, const uint target_nucleus_index,
                      const vector<spin>& nuclei, const double static_B, const int ms,
                      const uint k_DD, const double scale){
  // objects specific to target nucleus
  const spin target = nuclei.at(target_nucleus_index);
  const MatrixXcd rho_target_0 = MatrixXcd::Identity(2,2);
  const Vector3d target_A = A(target);
  const Vector3d target_A_perp = target_A - dot(target_A,zhat)*zhat;
  const Vector3d target_larmor = effective_larmor(target, static_B, ms);

  // minimum difference in larmor frequencies between target nucleus and other nuclei
  double dw_min = DBL_MAX;
  for(uint s = 0; s < nuclei.size(); s++){
    if(s == target_nucleus_index) continue;
    const double dw = abs(target_larmor.norm()
                          - effective_larmor(nuclei.at(s), static_B, ms).norm());
    if(dw < dw_min) dw_min = dw;
  }

  // AXY sequence parameters
  const double f_DD = 2*dw_min/(target_A_perp.norm()*scale);
  const double w_DD = target_larmor.norm()/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period

  // phase of hyperfine field rotation in the frame of H_nZ
  const double target_phi = acos(dot(hat(target_A_perp),xhat));

  // AXY sequence offsets for resonance with xhat and yhat componenets of the hyperfine field
  const double sx_pulse_offset = target_phi/(2*pi);
  const double sy_pulse_offset = (target_phi-0.5*pi)/(2*pi);

  // AXY protocol pulse times (normalized to one AXY sequence period)
  const vector<double> sx_pulses = axy_pulses(k_DD, f_DD, sx_pulse_offset);
  const vector<double> sy_pulses = axy_pulses(k_DD, f_DD, sy_pulse_offset);

  // NV center pulses
  const MatrixXcd X = act(sx,{0},2);
  const MatrixXcd Z_to_X = act(Ry(pi/2),{0},2);
  const MatrixXcd X_to_Z = act(Ry(-pi/2),{0},2);
  const MatrixXcd Z_to_Y = act(Rx(-pi/2),{0},2);
  const MatrixXcd Y_to_Z = act(Rx(pi/2),{0},2);

  // full NV+cluster Hamiltonian
  const MatrixXcd H = H_int_large_static_B(e(ms), {target}) + H_nZ({target}, static_B*zhat);

  MatrixXcd U_sx = exp(-j*H*t_DD*sx_pulses.at(0));
  for(uint i = 1; i < sx_pulses.size(); i++){
    U_sx = (exp(-j*H*t_DD*(sx_pulses.at(i)-sx_pulses.at(i-1))) * X * U_sx).eval();
  }
  U_sx = (exp(-j*H*t_DD*(1.-sx_pulses.back())) * X * U_sx).eval();
  U_sx = (Z_to_X * U_sx * X_to_Z).eval();

  MatrixXcd U_sy = exp(-j*H*t_DD*sy_pulses.at(0));
  for(uint i = 1; i < sy_pulses.size(); i++){
    U_sy = (exp(-j*H*t_DD*(sy_pulses.at(i)-sy_pulses.at(i-1))) * X * U_sy).eval();
  }
  U_sy = (exp(-j*H*t_DD*(1.-sy_pulses.back())) * X * U_sy).eval();
  U_sy = (Z_to_Y * U_sy * Y_to_Z).eval();



  U_sx = pow(U_sx,50);

  MatrixXcd H_sx = j*log(U_sx);
  remove_artifacts(H_sx);

  remove_artifacts(U_sx,1e-4);
  cout << U_sx << endl << endl;

  MatrixXcd H_int_eff = 0.25*ms*f_DD*target_A_perp.norm()*cos(target_phi)*tp(sx,sx);
  MatrixXcd H_nZ_eff = -act(dot(target.g*static_B*zhat - 0.5*ms*target_A, target.S) ,{1},2);
  MatrixXcd U_eff = exp(-j*(H_int_eff+H_nZ_eff)*t_DD);

  U_eff = pow(U_eff,50);
  remove_artifacts(U_eff,1e-4);

  cout << U_eff << endl << endl;


  return 1;
}
