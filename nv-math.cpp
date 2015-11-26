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

// group spins into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<uint>> get_index_clusters(const vector<spin>& spins,
                                        const double min_coupling_strength){

  vector<vector<uint>> ind_clusters; // all clusters
  vector<uint> clustered; // indices of spins we have already clustered

  // loop over indices of all spins
  for(uint i = 0; i < spins.size(); i++){
    if(!in_vector(i,clustered)){

      // initialize current cluster
      vector<uint> cluster;

      cluster.push_back(i);
      clustered.push_back(i);

      // loop over all indices of cluster
      for(uint ci = 0; ci < cluster.size(); ci++) {

        // loop over all spins indices greater than cluster(ci)
        for(uint k = cluster.at(ci)+1; k < spins.size(); k++){

          // if cluster(ci) and spins(k) are interacting, add k to this cluster
          if(!in_vector(k,clustered) &&
             coupling_strength(spins.at(cluster.at(ci)),
                               spins.at(k)) >= min_coupling_strength){
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
double find_target_coupling(const vector<spin>& spins, const uint cluster_size_target,
                            const double initial_cluster_coupling, const double dcc_cutoff){

  double cluster_coupling = initial_cluster_coupling;
  double dcc = cluster_coupling/4;

  vector<vector<uint>> ind_clusters = get_index_clusters(spins,cluster_coupling);
  bool coupling_too_small = largest_cluster_size(ind_clusters) >= cluster_size_target;
  bool last_coupling_too_small;
  bool crossed_correct_coupling;

  // determine coupling for which the largest cluster size just barely >= max_cluster_size
  while(dcc >= dcc_cutoff || !coupling_too_small){
    last_coupling_too_small = coupling_too_small;

    cluster_coupling += coupling_too_small ? dcc : -dcc;
    ind_clusters = get_index_clusters(spins,cluster_coupling);
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

// group spins into clusters according to cluster_indices
vector<vector<spin>> group_spins(const vector<spin>& spins,
                                 const vector<vector<uint>>& cluster_indices){
  vector<vector<spin>> clusters;
  for(uint c = 0; c < cluster_indices.size(); c++){
    vector<spin> cluster;
    for(uint s = 0; s < cluster_indices.at(c).size(); s++){
      cluster.push_back(spins.at(cluster_indices.at(c).at(s)));
    }
    clusters.push_back(cluster);
  }
  return clusters;
}

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// hyperfine field experienced by spin s
inline Vector3d A(const spin& s){
  const Vector3d r = s.pos - e(1).pos;
  return e(1).g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*dot(hat(r),zhat)*hat(r));
}

// effective larmor frequency of spin s
double effective_larmor(const spin& s, const Vector3d& B, const Vector3d& A, int ms){
  return (s.g*B - ms/2.*A).norm();
}

// pulse times for harmonic h and fourier component f
vector<double> pulse_spacings(const uint k, const double f){
  assert((k == 1) || (k == 3));
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

  vector<double> times;
  times.push_back(x1);
  times.push_back(x2);
  return times;
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
                             const double scan_time, const Vector3d& B){
  const spin e_ms = e(ms); // electron spin

  const double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period

  const vector<double> xs = pulse_spacings(k_DD,f_DD); // AXY protocol pulse times
  const double t1 = xs.at(0)*t_DD;
  const double t2 = xs.at(1)*t_DD;
  const double t3 = t_DD/4;

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    // projections onto |ms> and |0> NV states
    const MatrixXcd proj_m = act(up*up.adjoint(), {0}, clusters.at(c).size()+1); // |ms><ms|
    const MatrixXcd proj_0 = act(dn*dn.adjoint(), {0}, clusters.at(c).size()+1); // |0><0|

    // construct full Hamiltonian
    const MatrixXcd H = H_int(e_ms, clusters.at(c)) + H_Z(e_ms, clusters.at(c), B);

      // spin cluster Hamiltonians for single NV states
    const MatrixXcd H_m = ptrace(H*proj_m, {0}); // <ms|H|ms>
    const MatrixXcd H_0 = ptrace(H*proj_0, {0}); // <0|H|0>

    // propagators for sections of the AXY sequence
    const MatrixXcd U1_m = exp(-j*H_m*(t1));
    const MatrixXcd U2_m = exp(-j*H_0*(t2-t1));
    const MatrixXcd U3_m = exp(-j*H_m*(t3-t2));

    const MatrixXcd U1_0 = exp(-j*H_0*(t1));
    const MatrixXcd U2_0 = exp(-j*H_m*(t2-t1));
    const MatrixXcd U3_0 = exp(-j*H_0*(t3-t2));

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
                             double scan_time, const Vector3d& B_static,
                             const control_fields& controls, const uint integration_factor){
  const spin e_ms = e(ms); // electron spin

  const double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period

  double max_possible_B_mag = B_static.norm();
  for(uint c = 0; c < controls.num(); c++){
    max_possible_B_mag += controls.Bs.at(c).norm();
  }
  // todo: compute freq_scale, dt, and integration_steps for each cluster
  //   for the event that we simulate heteronuclear systems
  const double freq_scale = max(w_DD,gC13*max_possible_B_mag);
  const double dt = 1/(freq_scale*integration_factor); // integration time step

  scan_time = int(scan_time/t_DD)*t_DD; // scan an integer number of AXY sequences
  const uint integration_steps = int(scan_time/dt);

  const vector<double> xs = pulse_spacings(k_DD,f_DD); // AXY protocol pulse times
  const double t1 = xs.at(0)*t_DD;
  const double t2 = xs.at(1)*t_DD;
  const double t3 = t_DD/4;
  const double t4 = t_DD/2-t2;
  const double t5 = t_DD/2-t1;

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

    for(uint t_i = 1; t_i <= integration_steps; t_i++){
      const double t = t_i*dt; // time

      // flip electron spin according to the AXY sequence
      const double t_hAXY = t - floor(t/(t_DD/2))*t_DD/2; // time in current AXY half-sequence
      // if we are within dt/2 of an AXY pulse time, flip the projections
      if(min({abs(t_hAXY-t1), abs(t_hAXY-t2), abs(t_hAXY-t3),
              abs(t_hAXY-t4), abs(t_hAXY-t5)}) < dt/2){
        const MatrixXcd proj_tmp = proj_m;
        proj_m = proj_0;
        proj_0 = proj_tmp;
      }

      // compute net magnetic field
      Vector3d B = B_static;
      for(uint c = 0; c < controls.num(); c++){
        B += controls.Bs.at(c) * cos(controls.freqs.at(c)*t + controls.phases.at(c));
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

