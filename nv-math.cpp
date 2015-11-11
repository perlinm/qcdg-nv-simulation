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
double coupling_strength(const spin s1, const spin s2){
  Vector3d r = s2.pos - s1.pos;
  return s1.g*s2.g/(8*pi*pow(r.norm()*a0,3)) * (1 - 3*pow(hat(r).dot(zhat),2));
}

// group spins into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<spin>> get_clusters(vector<spin> spins, double min_coupling_strength){

  vector<vector<spin>> clusters; // vector of all spin clusters
  vector<uint> clustered; // indices of spins we have already clustered

  // loop over indices of all spins
  for(uint i = 0; i < spins.size(); i++){
    if(!in_vector(i,clustered)){

      // initialize current cluster
      vector<uint> cluster_indices;
      vector<spin> cluster;

      cluster_indices.push_back(i);
      cluster.push_back(spins.at(i));
      clustered.push_back(i);

      // loop over all indices of cluster
      for(uint ci = 0; ci < cluster_indices.size(); ci++) {

        // loop over all spins indices greater than cluster_indices(ci)
        for(uint k = cluster_indices.at(ci)+1; k < spins.size(); k++){

          // if cluster(ci) and spins(k) are interacting, add k to this cluster
          if(!in_vector(k,clustered) &&
             coupling_strength(cluster.at(ci),spins.at(k)) >= min_coupling_strength){
            cluster_indices.push_back(k);
            cluster.push_back(spins.at(k));
            clustered.push_back(k);
          }
        }
      }
      clusters.push_back(cluster);
    }
  }
  return clusters;
}

// get size of largest spin cluster
uint largest_cluster_size(vector<vector<spin>> clusters){
  uint largest_size = 0;
  for(uint i = 0; i < clusters.size(); i++){
    if(clusters.at(i).size() > largest_size) largest_size = clusters.at(i).size();
  }
  return largest_size;
}

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(vector<spin> spins, uint cluster_size_target,
                            double initial_cluster_coupling, double dcc_cutoff){

  double cluster_coupling = initial_cluster_coupling;
  double dcc = cluster_coupling/4;

  vector<vector<spin>> clusters = get_clusters(spins,cluster_coupling);
  bool coupling_too_small = largest_cluster_size(clusters) >= cluster_size_target;
  bool last_coupling_too_small;
  bool crossed_correct_coupling;

  // determine coupling for which the largest cluster size just barely >= max_cluster_size
  while(dcc >= dcc_cutoff || !coupling_too_small){
    last_coupling_too_small = coupling_too_small;

    cluster_coupling += coupling_too_small ? dcc : -dcc;
    clusters = get_clusters(spins,cluster_coupling);
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

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// hyperfine field experienced by spin s
Vector3d A(const spin& s){
  Vector3d r = s.pos - e(1).pos;
  return e(1).g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*dot(hat(r),zhat)*hat(r));
}

// effective larmor frequency of spin s
double effective_larmor(const spin& s, const Vector3d& B, const Vector3d& A, int ms){
  return (s.g*B - ms/2.*A).norm();
}

// pulse times for harmonic h and fourier component f
vector<double> pulse_spacings(uint k, double f){
  assert((k == 1) || (k == 3));
  double fp = f*pi;
  double x1,x2;
  if(k == 1){
    assert(abs(fp) < 8*cos(pi/9)-4);
    double w1 = 4 - fp;
    double w2 = w1 * (960 - 144*fp - 12*fp*fp + fp*fp*fp);

    x1 = 1/(2*pi) * atan2( (3*fp-12)*w1 + sqrt(3*w2),
                           sqrt(6)*sqrt(w2 - 96*fp*w1 + w1*w1*sqrt(3*w2)));
    x2 = 1/(2*pi) * atan2(-(3*fp-12)*w1 + sqrt(3*w2),
                          sqrt(6)*sqrt(w2 - 96*fp*w1 - w1*w1*sqrt(3*w2)));
  } else{ // if k == 3
    assert(abs(fp) < 4);
    double q1 = 4/(sqrt(5+fp)-1);
    double q2 = 4/(sqrt(5+fp)+1);

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
  Vector3d r = s2.pos-s1.pos;
  return s1.g*s2.g/(4*pi*pow(r.norm()*a0,3))
    * (dot(s1.S,s2.S) - 3*tp(dot(s1.S,hat(r)), (dot(s2.S,hat(r)))));
}

// return spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const spin& e, const vector<spin>& cluster){
  int spins = cluster.size()+1;
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

// return Zeeman Hamiltonian for NV center with cluster
MatrixXcd H_Z(const spin& e, const vector<spin>& cluster, const Vector3d& B){
  int spins = cluster.size()+1;
  // zero-field splitting and interaction of NV center with magnetic field
  MatrixXcd H = act(NV_ZFS*dot(e.S,zhat)*dot(e.S,zhat) - e.g*dot(B,e.S), {0}, spins);
  for(uint s = 0; s < cluster.size(); s++){
    // interaction of spin s with magnetic field
    H -= act(cluster.at(s).g*dot(B,cluster.at(s).S), {s+1}, spins);
  }
  return H;
}

// perform NV coherence measurement
double coherence_measurement(int ms, const vector<vector<spin>>& clusters, double w_scan,
                             uint k_DD, double f_DD, const Vector3d& B, double scan_time){
  spin e_ms = e(ms); // electron spin

  double w_DD = w_scan/k_DD; // AXY protocol angular frequency
  double t_DD = 2*pi/w_DD; // AXY protocol period

  vector<double> xs = pulse_spacings(k_DD,f_DD); // AXY protocol pulse times
  double t1 = xs.at(0)*t_DD;
  double t2 = xs.at(1)*t_DD;
  double t3 = t_DD/4;

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    // projections onto |ms> and |0> NV states
    MatrixXcd proj_m = act(up*up.adjoint(), {0}, clusters.at(c).size()+1); // |ms><ms|
    MatrixXcd proj_0 = act(dn*dn.adjoint(), {0}, clusters.at(c).size()+1); // |0><0|

    // construct full Hamiltonian
    MatrixXcd H = H_int(e_ms, clusters.at(c)) + H_Z(e_ms, clusters.at(c), B);

    // extract sub-system Hamiltonians
    MatrixXcd H_m = ptrace(H*proj_m, {0});
    MatrixXcd H_0 = ptrace(H*proj_0, {0});

    // propagators for sections of the AXY sequence
    MatrixXcd U1_m = exp(-j*H_m*(t1));
    MatrixXcd U2_m = exp(-j*H_0*(t2-t1));
    MatrixXcd U3_m = exp(-j*H_m*(t3-t2));

    MatrixXcd U1_0 = exp(-j*H_0*(t1));
    MatrixXcd U2_0 = exp(-j*H_m*(t2-t1));
    MatrixXcd U3_0 = exp(-j*H_0*(t3-t2));

    MatrixXcd U4_m = U3_0;
    MatrixXcd U5_m = U2_0;
    MatrixXcd U6_m = U1_0;

    MatrixXcd U4_0 = U3_m;
    MatrixXcd U5_0 = U2_m;
    MatrixXcd U6_0 = U1_m;

    // single AXY sequence propagators
    MatrixXcd U_m = U1_m*U2_m*U3_m*U4_m*U5_m*U6_m * U6_m*U5_m*U4_m*U3_m*U2_m*U1_m;
    MatrixXcd U_0 = U1_0*U2_0*U3_0*U4_0*U5_0*U6_0 * U6_0*U5_0*U4_0*U3_0*U2_0*U1_0;

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
//--------------------------------------------------------------------------------------------

