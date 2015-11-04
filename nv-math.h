#pragma once

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "qp-math.h"
#include "constants.h"

// spin up/down state vectors
const VectorXcd up = (Vector2cd() << 1,0).finished();
const VectorXcd dn = (Vector2cd() << 0,1).finished();

// pauli spin matrices
const MatrixXcd st = (Matrix2cd() << 1,0, 0,1).finished();
const MatrixXcd sx = (Matrix2cd() << 0,1, 1,0).finished();
const MatrixXcd sy = (Matrix2cd() << 0,-j, j,0).finished();
const MatrixXcd sz = (Matrix2cd() << 1,0, 0,-1).finished();

// spin vector for a spin-1/2 particle
const mvec s_vec({ sx/2, sy/2, sz/2 });

// diamond lattice vectors scaled to a0
const Vector3d ao = (Vector3d() << 1,1,1).finished()/4;
const Vector3d a1 = (Vector3d() << 0,1,1).finished()/2;
const Vector3d a2 = (Vector3d() << 1,0,1).finished()/2;
const Vector3d a3 = (Vector3d() << 1,1,0).finished()/2;

// vector of lattice sites in a diamond unit cell
const vector<Vector3d> cell_sites { Vector3d::Zero(), a1, a2, a3, ao, ao+a1, ao+a2, ao+a3 };

// struct for spins
struct spin{
  Vector3d pos; // position
  double g; // gyromagnetic ratio
  mvec S; // spin vector

  spin(Vector3d pos, double g, vector<MatrixXcd> S){
    assert(S.size() == 3);
    this->pos = pos;
    this->g = g;
    this->S = mvec(S);
  };
  spin(Vector3d pos, double g, mvec S){
    assert(S.size() == 3);
    this->pos = pos;
    this->g = g;
    this->S = S;
  };

  bool operator==(spin s){
    return ((pos == s.pos) && (g == s.g) && (S == s.S));
  }
  bool operator!=(spin s){ return !(*this == s); }

};

// initialize nitrogen and vacancy centers
const spin n(Vector3d::Zero(), 0., s_vec);
spin e(int ms){
  return spin(ao, ge, {sx/sqrt(2), ms*sy/sqrt(2), ms*(sz+I2)/2.});
}

// unit vectors along lattice directions
const Vector3d zhat = (Vector3d() << 1,1,1).finished()/sqrt(3); // direction from N to V site
const Vector3d xhat = (Vector3d() << 2,-1,-1).finished()/sqrt(6);
const Vector3d yhat = (Vector3d() << 0,1,-1).finished()/sqrt(2);

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// coupling strength between two spins; assumes strong magnetic field in zhat
inline double coupling_strength(const spin s1, const spin s2){
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

// hyperfine field experienced by nuclear spin s
Vector3d A(const spin s, int ms){
  Vector3d r = s.pos - e(ms).pos;
  return e(ms).g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*hat(r).dot(zhat)*hat(r));
}

// harmonic to target with AXY sequence
enum harmonic { first, third };

// pulse times for harmonic h and fourier component f
vector<double> pulse_times(harmonic k, double f){

  double fp = f*pi;
  double t1,t2;
  if(k == first){
    assert(abs(fp) < 8*cos(pi/9)-4);
    double w1 = 4 - fp;
    double w2 = w1 * (960 - 144*fp - 12*fp*fp + fp*fp*fp);

    t1 = 1/(2*pi) * atan2( (3*fp-12)*w1 + sqrt(3*w2),
                           sqrt(6)*sqrt(w2 - 96*fp*w1 + w1*w1*sqrt(3*w2)));
    t2 = 1/(2*pi) * atan2(-(3*fp-12)*w1 + sqrt(3*w2),
                          sqrt(6)*sqrt(w2 - 96*fp*w1 - w1*w1*sqrt(3*w2)));
  } else{ // if harmonic == third
    assert(abs(fp) < 4);
    double q1 = 4/(sqrt(5+fp)-1);
    double q2 = 4/(sqrt(5+fp)+1);

    t1 = 1/4 - 1/(2*pi)*atan(sqrt(q1*q1-1));
    t2 = 1/4 - 1/(2*pi)*atan(sqrt(q2*q2-1));
  }

  vector<double> times;
  times.push_back(t1);
  times.push_back(t2);
  return times;
}

// Hamiltoninan coupling two spins
inline MatrixXcd H_ss(const spin s1, const spin s2){
  Vector3d r = s2.pos-s1.pos;
  return s1.g*s2.g/(4*pi*pow(r.norm()*a0,3))* (dot(s1.S,s2.S)
                                               - 3*tp(dot(s1.S,hat(r)), (dot(s2.S,hat(r)))));
}

// computer NV coherence at scanning frequency w_scan
double coherence_scan(vector<vector<spin>> clusters, double w_scan, harmonic k_DD,
                      double f_DD, double Bz, int ms, double scan_time){
  double w_DD = w_scan;
  if(k_DD == third) w_DD /= 3;

  double t_DD = 2*pi/w_DD;

  vector<double> ts = pulse_times(k_DD,f_DD);
  double t1 = ts.at(0)*t_DD;
  double t2 = ts.at(1)*t_DD;
  double t3 = t_DD/4;

  // initial (unnormalized) density matrix of NV center
  // todo: do this in the NV z-axis basis
  MatrixXcd rho_NV_0 = (up+dn)*(up+dn).adjoint();

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    vector<spin> cluster = clusters.at(c);
    int spins = cluster.size()+1; // total number of spins in NV+cluster system
    int D = pow(2,spins); // dimensionality of NV+cluster Hilbert space

    // initial (unnormalized) density matrix of NV+cluster
    MatrixXcd rho_0 = act(rho_NV_0,{0},spins);

    // pi-pulse on NV center
    MatrixXcd X = act(sx,{0},spins);

    // construct Hamiltonians
    MatrixXcd H_nn = MatrixXcd::Zero(D,D); // internuclear coupling Hamiltonian
    MatrixXcd H_nZ_eff = MatrixXcd::Zero(D,D); // effective nuclear Zeeman Hamiltonian
    MatrixXcd H_int = MatrixXcd::Zero(D,D); // NV-cluster interaction Hamiltonian

    // loop over spins in cluster
    for(uint s = 0; s < cluster.size(); s++){
      for(uint r = 0; r < s; r++){
        H_nn += act(H_ss(cluster.at(r), cluster.at(s)), {r+1,s+1}, spins);
      }
      MatrixXcd AI = dot(A(cluster.at(s), ms), cluster.at(s).S);
      H_nZ_eff += act( -cluster.at(s).g*dot(Bz*zhat, cluster.at(s).S) + ms/2.*AI,
                       {s+1}, spins);
      /* H_int += ms/2. * act(tp(dot(e(ms).S,zhat), AI), {0,s+1}, spins); */
      H_int += ms/2. * act(tp(sz, AI), {0,s+1}, spins);
    }

    MatrixXcd H = H_nn + H_nZ_eff + H_int;

    // propagators for sections of the AXY sequence
    MatrixXcd U1 = exp(-j*H*(t1));
    MatrixXcd U2 = exp(-j*H*(t2-t1));
    MatrixXcd U3 = exp(-j*H*(t3-t2));

    // AXY half-sequence propagator
    MatrixXcd U = U1*X*U2*X*U3*X*U3*X*U2*X*U1;

    // propagator for entire scan
    U = pow(U,2*int(scan_time/t_DD));

    // (unnormalized) NV+cluster density matrix after scanning
    MatrixXcd rho = U*rho_0*U.adjoint();

    // update coherence
    coherence *= real(trace(rho*rho_0))/D-1;
  }
  return coherence;
}

// return spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const spin e, const vector<spin> cluster){
  int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));

  for(uint s = 0; s < cluster.size(); s++){
    // interaction with NV center
    H += act(H_ss(e, cluster.at(s)), {0,s+1}, spins);
    for(uint r = 0; r < s; r++){
      // interaction with other spins
      H += act(H_ss(cluster.at(r), cluster.at(s)), {r+1,s+1}, spins);
    }
  }
  return H;
}

// return Zeeman Hamiltonian for NV center with cluster
MatrixXcd H_Z(const spin e, const vector<spin> cluster, Vector3d B){
  int spins = cluster.size()+1;

  MatrixXcd H = act(NV_ZFS*dot(e.S,zhat)*dot(e.S,zhat) - e.g*dot(B,e.S), {0}, spins);
  for(uint s = 0; s < cluster.size(); s++){
    H -= act(cluster.at(s).g*dot(B,cluster.at(s).S), {s+1}, spins);
  }
  return H;
}

// computer NV coherence at scanning frequency w_scan
double exact_coherence_scan(vector<vector<spin>> clusters, double w_scan, harmonic k_DD,
                            double f_DD, double Bz, int ms, double scan_time){
  spin e_ms = e(ms);
  Vector3d B = Bz*zhat;

  double w_DD = w_scan;
  if(k_DD == third) w_DD /= 3;

  double t_DD = 2*pi/w_DD;

  vector<double> ts = pulse_times(k_DD,f_DD);
  double t1 = ts.at(0)*t_DD;
  double t2 = ts.at(1)*t_DD;
  double t3 = t_DD/4;

  // initial (unnormalized) density matrix of NV center
  MatrixXcd rho_NV_0 = (up+dn)*(up+dn).adjoint();

  double coherence = 1;
  for(uint c = 0; c < clusters.size(); c++){
    vector<spin> cluster = clusters.at(c);
    int spins = cluster.size()+1;

    // pi-pulse on NV center
    MatrixXcd X = act(sx,{0},spins);

    // initial (unnormalized) density matrix of NV+cluster
    MatrixXcd rho_0 = act(rho_NV_0,{0},spins);

    // construct full Hamiltonian
    MatrixXcd H = H_int(e_ms, cluster) + H_Z(e_ms, cluster, B);

    // propagators for sections of the AXY sequence
    MatrixXcd U1 = exp(-j*H*(t1));
    MatrixXcd U2 = exp(-j*H*(t2-t1));
    MatrixXcd U3 = exp(-j*H*(t3-t2));

    // AXY half-sequence propagator
    MatrixXcd U = U1*X*U2*X*U3*X*U3*X*U2*X*U1;

    // propagator for entire scan
    U = pow(U,2*int(scan_time/t_DD));

    // (unnormalized) NV+cluster density matrix after scanning
    MatrixXcd rho = U*rho_0*U.adjoint();

    // update coherence
    coherence *= real(trace(rho*rho_0))/pow(2,spins)-1;
  }
  return coherence;
}
