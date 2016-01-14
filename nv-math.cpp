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

nv_system::nv_system(const int ms, const uint k_DD,
                     const double static_Bz, const double scale_factor) :
  e(spin(Vector3d::Zero(), ge,
         mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat))),
  ms(ms), k_DD(k_DD), static_Bz(static_Bz), scale_factor(scale_factor)
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
vector<vector<uint>> get_index_clusters(const vector<spin>& nuclei,
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

  vector<vector<uint>> clusters = get_index_clusters(nuclei, cluster_coupling);
  bool coupling_too_small = largest_cluster_size(clusters) >= cluster_size_target;
  bool last_coupling_too_small;
  bool crossed_correct_coupling = false;

  // determine coupling for which the largest cluster size just barely >= max_cluster_size
  while(dcc >= dcc_cutoff || !coupling_too_small){
    last_coupling_too_small = coupling_too_small;

    cluster_coupling += coupling_too_small ? dcc : -dcc;
    clusters = get_index_clusters(nuclei,cluster_coupling);
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

// convert cluster of indices into cluster of spins
vector<vector<spin>> spin_clusters(const vector<spin>& nuclei,
                                   const vector<vector<uint>>& index_clusters){
  vector<vector<spin>> spin_clusters;
  for(uint c = 0; c < index_clusters.size(); c++){
    vector<spin> spin_cluster;
    for(uint s = 0; s < index_clusters.at(c).size(); s++){
      spin_cluster.push_back(nuclei.at(index_clusters.at(c).at(s)));
    }
    spin_clusters.push_back(spin_cluster);
  }
  return spin_clusters;
}

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// pulse times for harmonic h and fourier component f
vector<double> axy_pulses(const uint k, const double f){

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
  return pulses;
}

vector<double> delayed_pulses(const vector<double> pulses, const double delay){
  // number of pulses
  const uint N = pulses.size();

  // delayed pulses
  vector<double> delayed_pulses;
  for(uint p = 0; p < 2*N; p++){
    if(p/N + pulses.at(p%N) - delay > 0){
      delayed_pulses.push_back(p/N + pulses.at(p%N) - delay);
    }
    if(delayed_pulses.size() == N) break;
  }
  return delayed_pulses;
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
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time){
  const vector<vector<spin>> clusters = spin_clusters(nv.nuclei,nv.clusters); // spin clusters

  const double w_DD = w_scan/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulses = axy_pulses(nv.k_DD,f_DD); // AXY protocol pulse times

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    // projections onto |ms> and |0> NV states
    const MatrixXcd proj_m = act(up*up.adjoint(), {0}, clusters.at(c).size()+1); // |ms><ms|
    const MatrixXcd proj_0 = act(dn*dn.adjoint(), {0}, clusters.at(c).size()+1); // |0><0|

    // construct full Hamiltonian
    const MatrixXcd H = H_int(nv.e,clusters.at(c)) + H_Z(nv.e,clusters.at(c),
                                                         nv.static_Bz*zhat);

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

    // normalize propagators
    U_m /= sqrt(real(trace(U_m.adjoint()*U_m)/double(U_m.rows())));
    U_0 /= sqrt(real(trace(U_0.adjoint()*U_0)/double(U_0.rows())));

    // update coherence
    coherence *= real(trace(U_0.adjoint()*U_m)) / pow(2,clusters.at(c).size());
  }
  return coherence;
}

//--------------------------------------------------------------------------------------------
// Control field scanning
// (assume large static magnetic field along NV axis)
//--------------------------------------------------------------------------------------------

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int_large_static_B(const nv_system& nv, const vector<spin>& cluster){
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for(uint s = 0; s < cluster.size(); s++){
    // interaction between NV center and spin s
    H += act( tp(dot(nv.e.S,zhat), dot(A(nv,cluster.at(s)),cluster.at(s).S)), {0,s+1}, spins);
    for(uint r = 0; r < s; r++){
      // interaction betwen spin r and spin s
      H += act(H_ss_large_static_B(cluster.at(r), cluster.at(s)), {r+1,s+1}, spins);
    }
  }
  return H;
}

// perform NV coherence measurement with a static magnetic field and additional control fields
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time, const control_fields& controls,
                             const uint integration_factor){
  const vector<vector<spin>> clusters = spin_clusters(nv.nuclei,nv.clusters); // spin clusters

  const double w_DD = w_scan/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulses = axy_pulses(nv.k_DD,f_DD); // AXY protocol pulse times

  Vector3d max_possible_B = nv.static_Bz*zhat;
  for(uint c = 0; c < controls.num(); c++){
    Vector3d B = controls.Bs.at(c);
    max_possible_B += abs(dot(B,xhat))*xhat + abs(dot(B,yhat))*yhat + abs(dot(B,zhat))*zhat;
  }

  // todo: compute freq_scale, dt, dx, and integration_steps for each cluster
  //   for the event that we simulate heteronuclear systems
  const double freq_scale = max(w_DD,gC13*max_possible_B.norm());
  const double dt = 1/(freq_scale*integration_factor);
  const double dx = 1/(t_DD*freq_scale*integration_factor);

  // use an integer number of integration steps to scan an integer number of AXY sequences
  const uint integration_steps = int(int(scan_time/t_DD)*t_DD / dt);

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

      // net magnetic field
      Vector3d B = nv.static_Bz*zhat;
      for(uint c = 0; c < controls.num(); c++){
        B += controls.Bs.at(c) * cos(controls.freqs.at(c)*x*t_DD + controls.phases.at(c));
      }

      // NV+cluster Hamiltonian
      const MatrixXcd H = H_int_large_static_B(nv,clusters.at(c)) + H_nZ(clusters.at(c),B);

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
fidelity_info iswap_fidelity(const nv_system& nv, const uint index){
  // identify cluster of target nucleus
  const vector<vector<spin>> clusters = spin_clusters(nv.nuclei,nv.clusters);
  vector<spin> cluster;
  uint cluster_index;
  for(uint c = 0; c < clusters.size(); c++){
    for(uint s = 0; s < clusters.at(c).size(); s++){
      if(clusters.at(c).at(s) == nv.nuclei.at(index)){
        cluster = clusters.at(c);
        cluster_index = s;
      }
    }
  }

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
  const double w_DD = larmor_eff.norm()/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const double operation_time = 2*pi/abs(f_DD*hyperfine_perp.norm());

  // AXY pulse sequence matching target larmor frequency
  const vector<double> sx_pulses = axy_pulses(nv.k_DD, f_DD);
  const vector<double> sy_pulses = delayed_pulses(axy_pulses(nv.k_DD, f_DD), 0.25);

  // NV+cluster Hamiltonian
  const MatrixXcd H = H_int(nv.e, cluster) + H_Z(nv.e, cluster, nv.static_Bz*zhat);

  // NV center spin flip pulse
  const MatrixXcd X = act(sx,{0},cluster.size()+1);

  // spin addressing propagators
  MatrixXcd U_sx = X * exp(-j*H*t_DD*sx_pulses.at(0));
  MatrixXcd U_sy = X * exp(-j*H*t_DD*sy_pulses.at(0));
  for(uint i = 1; i < sx_pulses.size(); i++){
    U_sx = (X * exp(-j*H*t_DD*(sx_pulses.at(i)-sx_pulses.at(i-1))) * U_sx).eval();
    U_sy = (X * exp(-j*H*t_DD*(sy_pulses.at(i)-sy_pulses.at(i-1))) * U_sy).eval();
  }
  U_sx = (exp(-j*H*t_DD*(1.-sx_pulses.back())) * U_sx).eval();
  U_sy = (exp(-j*H*t_DD*(1.-sy_pulses.back())) * U_sy).eval();

  U_sx = pow(U_sx, int(operation_time/t_DD));
  U_sy = pow(U_sy, int(operation_time/t_DD));

  // "natural" basis vectors for target nucleus
  const Vector3d target_zhat = hat(larmor_eff);
  const Vector3d target_xhat = hat(hyperfine_perp);
  const Vector3d target_yhat = target_zhat.cross(target_xhat);

  // operators to rotate into the "natural" frame
  const Matrix2cd Z_to_tY = rotate(acos(dot(zhat,target_yhat)), zhat.cross(target_yhat));
  const Matrix2cd Z_to_tX = rotate(acos(dot(zhat,target_xhat)), zhat.cross(target_xhat));

  // iSWAP operation
  const MatrixXcd iSWAP =
    act(Z_to_tX,{0},cluster.size()+1) * U_sx * act(Z_to_tX.inverse(),{0},cluster.size()+1) *
    act(Z_to_tY,{0},cluster.size()+1) * U_sy * act(Z_to_tY.inverse(),{0},cluster.size()+1);

  // spin operators in the "natural" basis
  const Matrix2cd sxp = dot(s_vec,target_xhat);
  const Matrix2cd syp = dot(s_vec,target_yhat);

  // exact iSWAP gate
  const MatrixXcd iSWAP_exact = exp(j*0.25*pi * act( tp(sxp,sxp) + tp(syp,syp),
                                                    {0,cluster_index+1}, cluster.size()+1));

  // compute iSWAP gate fidelity
  const double fidelity = gate_fidelity(iSWAP,iSWAP_exact);

  return fidelity_info(larmor_eff.norm(), hyperfine.norm(), hyperfine_perp.norm(),
                       dw_min, f_DD, operation_time, fidelity);
}

// realization of propagator U = exp(-i * phi * sigma_{axis}^{index})
MatrixXcd U_ctl(const nv_system& nv, const uint index, const double phi, const Vector3d axis,
                const uint interaction_factor, const double threshold){
  assert(abs(dot(hat(axis),zhat)) < threshold);

  const vector<vector<spin>> clusters = spin_clusters(nv.nuclei,nv.clusters);
  vector<spin> cluster;
  uint cluster_index;
  for(uint c = 0; c < clusters.size(); c++){
    for(uint s = 0; s < clusters.at(c).size(); s++){
      if(clusters.at(c).at(s) == nv.nuclei.at(index)){
        cluster = clusters.at(c);
        cluster_index = s;
      }
    }
  }
  const uint spins = cluster.size()+1;

  // control field frequency = effective larmor frequency of target nucleus
  const double larmor_eff = effective_larmor(nv,index).norm();

  double w_min = DBL_MAX; // minimum larmor frequency of all nuclei
  double dw_min = DBL_MAX; // min{ |w_s - w_{target}| for all s }
  for(uint s = 0; s < nv.nuclei.size(); s++){
    const double larmor_eff_s = effective_larmor(nv,s).norm();
    if(larmor_eff_s < w_min) w_min = larmor_eff_s;
    if(s == index) continue;
    const double dw = abs(larmor_eff - larmor_eff_s);
    if(dw < dw_min) dw_min = dw;
  }

  // AXY protocol frequency and period
  const double w_DD = w_min/2;
  const double t_DD = 2*pi/w_DD;

  // todo: try w_ctl = gC13*nv.static_Bz/2;
  const double w_ctl = larmor_eff; // control field frequency
  const double g_B = dw_min/nv.scale_factor; // control field strength * gyromagnetic ratio
  const double operation_time = 4*phi/g_B; // control operation time

  const double dt = 1/(w_ctl*integration_factor); // integration step size
  const uint integration_steps = int(operation_time/dt); // number of integration steps

  // initial propagator
  const MatrixXcd Id = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  MatrixXcd U = Id;

  for(uint t_i = 1; t_i < integration_steps; t_i++){
    const double t = t_i*dt;

    const MatrixXcd H = H_int_large_static_B(nv,cluster) + H_nZ(cluster,B);
    // flip electron spin according to the AXY sequence
    const double t_hAXY = t - floor(t/0.5)*0.5; // time in current AXY half-sequence
    // if we are within dx/2 of an AXY pulse time, flip the projections
    if(min({abs(t_hAXY-pulses.at(0)), abs(t_hAXY-pulses.at(1)), abs(t_hAXY-pulses.at(2)),
            abs(t_hAXY-pulses.at(3)), abs(t_hAXY-pulses.at(4))}) < dt/2){
      U *= act(X, {0}, spins);
    }

    // control field and Hamiltonian
    const Vector3d B = g_B/nv.nuclei.at(index).g * cos(w_ctl*t) * hat(axis);
    const MatrixXcd H = H_int_large_static_B(nv,clusters.at(c)) + H_nZ(clusters.at(c),B);

    // update and normalize propagator
    // U = (exp(-j*dt*H)*U).eval();
    // U = ((Id - j*dt*H - dt*dt*H*H/2)*U).eval();
    U = ((Id - j*dt*H)*U).eval();
    U /= sqrt(real(trace(U.adjoint()*U)/cHD));
  }

  return U;
}
