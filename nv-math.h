#pragma once

#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
#include <eigen3/unsupported/Eigen/KroneckerProduct> // provides tensor product
#include <eigen3/unsupported/Eigen/MatrixFunctions> // provides matrix functions
using namespace Eigen;

#include "constants.h"

// random double from 0 to 1
inline double rnd(){ return std::rand()/(double)RAND_MAX; }

// check whether value val is in vector vec
inline bool in_vector(auto val, vector<auto> vec){
  return (find(vec.begin(), vec.end(), val) != vec.end());
}

//--------------------------------------------------------------------------------------------
// Mathematical constants and methods
//--------------------------------------------------------------------------------------------

// identity matrices
const MatrixXcd I1 = MatrixXcd::Identity(1,1);
const MatrixXcd I2 = MatrixXcd::Identity(2,2);
const MatrixXcd I4 = MatrixXcd::Identity(4,4);

// return unit vector in direction of vec
inline Vector3d normed(Vector3d vec){ return vec/vec.norm(); }

// matrix functions
inline complex<double> trace(const MatrixXcd M){ return M.trace(); }
inline MatrixXcd log(const MatrixXcd M){ return M.log(); }
inline MatrixXcd exp(const MatrixXcd M){ return M.exp(); }
inline MatrixXcd sqrt(const MatrixXcd M){ return M.sqrt(); }
inline MatrixXcd pow(const MatrixXcd M, auto x){ return M.pow(x); }

// tensor product of two matrices
inline MatrixXcd tp(const MatrixXcd A, const MatrixXcd B){ return kroneckerProduct(A,B); }

// tensor product of many matrices
MatrixXcd tp(const initializer_list<MatrixXcd> list){
  MatrixXcd out = I1;
  for(MatrixXcd elem: list){
    out = tp(out,elem);
  }
  return out;
}

// inner product of two matrix vectors
inline MatrixXcd dot(const vector<MatrixXcd> v, const vector<MatrixXcd> w){
  assert(v.size() == w.size());
  MatrixXcd out = tp(v.at(0),w.at(0));
  for(int i = 1; i < v.size(); i++){
    out += tp(v.at(i),w.at(i));
  }
  return out;
}

// inner product of a matrix vector and a spacial vector
inline MatrixXcd dot(const vector<MatrixXcd> v, const Vector3d r){
  assert(v.size() == 3);
  return v.at(0)*r(0) + v.at(1)*r(1) + v.at(2)*r(2);
}
inline MatrixXcd dot(const Vector3d r, const vector<MatrixXcd> v){ return dot(v,r); }

// remove numerical artifacts from a matrix
void remove_artifacts(MatrixXcd &A, double threshold = 1e-12){
  for(int m = 0; m < A.rows(); m++){
    for(int n = 0; n < A.cols(); n++){
      if(abs(A(m,n).real()) < threshold) A(m,n) -= A(m,n).real();
      if(abs(A(m,n).imag()) < threshold) A(m,n) -= A(m,n).imag()*j;
    }
  }
}

// get global phase of matrix
complex<double> get_phase(const MatrixXcd A){
  for(int m = 0; m < A.rows(); m++){
    for(int n = 0; n < A.cols(); n++){
      if(abs(A(m,n)) != 0){
        complex<double> phase = A(m,n)/abs(A(m,n));
        if(phase.real() < 0) phase = -phase;
        return phase;
      }
    }
  }
  return 1;
}

// remove global phase from matrix
void remove_phase(MatrixXcd &A){
  A *= conj(get_phase(A));
}

//--------------------------------------------------------------------------------------------
// Spin matrix and propagator methods
//--------------------------------------------------------------------------------------------

// spin up/down state vectors
const VectorXcd up = (Vector2cd() << 1,0).finished();
const VectorXcd dn = (Vector2cd() << 0,1).finished();

// pauli spin matrices
const MatrixXcd st = (Matrix2cd() << 1,0, 0,1).finished();
const MatrixXcd sx = (Matrix2cd() << 0,1, 1,0).finished();
const MatrixXcd sy = (Matrix2cd() << 0,-j, j,0).finished();
const MatrixXcd sz = (Matrix2cd() << 1,0, 0,-1).finished();

// spin vector for a spin-1/2 particle
const vector<MatrixXcd> s_vec = { sx/2, sy/2, sz/2 };

// get the n-th bit of an integer num
inline bool int_bit(int num, int n){
  if(pow(2,n) > num) return 0;
  else return (num >> n) & 1;
}

// get state of qbit p (of N) from enumerated state s
inline bool qbit_state(int p, int N, int s){ return int_bit(s,N-1-p); }

// generate matrix B to act A on qbits qs_act out of qbits_new
MatrixXcd act(const MatrixXcd A, const vector<int> qs_act, int qbits_new){
  assert(A.rows() == A.cols()); // A should be square

  // number of qbits A acted on
  const int qbits_old = qs_act.size();
  assert(qbits_old == log2(A.rows()));

  // vector of qubits we are ignoring
  vector<int> qs_ignore;
  for(int i = 0; i < qbits_new; i++){
    if(!in_vector(i,qs_act)){
      qs_ignore.push_back(i);
    }
  }

  // initialize B (output) to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries A(m,n)
  for(int m = 0; m < A.rows(); m++){
    for(int n = 0; n < A.cols(); n++){

      // contribution of m and n to indices of B
      int b_m = 0, b_n = 0;
      for(int q = 0; q < qbits_old; q++){
        if(qbit_state(q,qbits_old,m)) b_m += pow(2,qbits_new-1-qs_act.at(q));
        if(qbit_state(q,qbits_old,n)) b_n += pow(2,qbits_new-1-qs_act.at(q));
      }

      // loop over all elements of the form |ms><ns| in B
      for(int s = 0; s < pow(2,qs_ignore.size()); s++){
        int b_out = b_m, b_in = b_n;
        for(int q = 0; q < qs_ignore.size(); q++){
          if(int_bit(s,q)){
            b_out += pow(2,qbits_new-1-qs_ignore.at(q));
            b_in += pow(2,qbits_new-1-qs_ignore.at(q));
          }
        }
        B(b_out,b_in) = A(m,n);
      }
    }
  }
  return B;
}

// perform a partial trace over qbits qs_trace
MatrixXcd ptrace(const MatrixXcd A, const vector<int> qs_trace){
  assert(A.rows() == A.cols()); // A should be square

  // number of qbits A acted on
  const int qbits_old = log2(A.rows());
  const int qbits_new = qbits_old - qs_trace.size();

  // vector of qubits we are keeping
  vector<int> qs_keep;
  for(int i = 0; i < qbits_old; i++){
    if(!in_vector(i,qs_trace)) qs_keep.push_back(i);
  }
  assert(qbits_new == qs_keep.size());

  // initialize B (output) to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries B(m,n)
  for(int m = 0; m < B.rows(); m++){
    for(int n = 0; n < B.cols(); n++){

      // contribution of m and n to indices of A
      int a_m = 0, a_n = 0;
      for(int q = 0; q < qbits_new; q++){
        if(qbit_state(q,qbits_new,m)) a_m += pow(2,qbits_old-1-qs_keep.at(q));
        if(qbit_state(q,qbits_new,n)) a_n += pow(2,qbits_old-1-qs_keep.at(q));
      }

      // loop over all elements of the form |ms><ns| in A
      for(int s = 0; s < pow(2,qs_trace.size()); s++){
        int a_out = a_m, a_in = a_n;
        for(int q = 0; q < qs_trace.size(); q++){
          if(int_bit(s,q)){
            a_out += pow(2,qbits_old-1-qs_trace.at(q));
            a_in += pow(2,qbits_old-1-qs_trace.at(q));
          }
        }
        B(m,n) += A(a_out,a_in);
      }
    }
  }
  return B;
}

//--------------------------------------------------------------------------------------------
// Qbit and matrix printing methods
//--------------------------------------------------------------------------------------------

// clean up matrix for human readability
inline void clean(MatrixXcd &M){
  remove_artifacts(M);
  remove_phase(M);
  remove_artifacts(M);
}

// returns basis element p for Hamiltonians of a system with N spins
MatrixXcd H_basis_element(int p, int N){
  MatrixXcd spins[4] = {st, sx, sy, sz};
  MatrixXcd b_p = spins[int_bit(p,0)+2*int_bit(p,1)];
  for(int n = 1; n < N; n++){
    b_p = tp(b_p,spins[int_bit(p,2*n)+2*int_bit(p,2*n+1)]);
  }
  return b_p;
}
string H_basis_element_text(int p, int N){
  char spins[4] = {'I','X','Y','Z'};
  stringstream stream;
  for(int n = 0; n < N; n++){
    stream << spins[int_bit(p,2*n)+2*int_bit(p,2*n+1)];
  }
  return stream.str();
}

// flatten matrix into a 1-D vector
inline MatrixXcd flatten(MatrixXcd M){
  M.resize(M.size(),1);
  return M;
}

// returns matrix whose columns are basis Hamiltonians for a system of N spins
MatrixXcd H_basis_matrix(int N){
  MatrixXcd spins[4] = {st, sx, sy, sz};
  MatrixXcd out = MatrixXcd::Zero(pow(4,N),pow(4,N));
  for(int p = 0; p < pow(4,N); p++){
    out.col(p) = flatten(H_basis_element(p,N));
  }
  return out;
}

// decompose Hamiltonian into its basis elements
MatrixXcd H_decompose(const MatrixXcd H, bool fast = true){
  int N = log2(H.rows());
  if(fast) return H_basis_matrix(N).householderQr().solve(flatten(H));
  else return H_basis_matrix(N).fullPivLu().solve(flatten(H));
}

// print Hamiltonian in human-readable form
void H_print(const MatrixXcd H, bool fast = true){
  int N = log2(H.rows());
  MatrixXcd hs = H_decompose(H,fast);
  clean(hs);
  for(int p = 0; p < pow(4,N); p++){
    if(abs(hs(p)) != 0){
      cout << H_basis_element_text(p,N) << ": " << hs(p) << endl;
    }
  }
}

// print state vector in human readable form
void state_print(const MatrixXcd psi){
  int N = psi.size();
  int qbits = log2(N);
  for(int n = 0; n < N; n++){
    if(abs(psi(n)) != 0){
      cout << "|";
      for(int q = 0; q < qbits; q++){
        cout << (qbit_state(n,qbits,q)?"d":"u");
      }
      cout << "> " << psi(n) << endl;
    }
  }
}

// print matrix in human readable form
void matrix_print(const MatrixXcd M){
  int qbits = log2(M.rows());
  for(int m = 0; m < M.rows(); m++){
    for(int n = 0; n < M.cols(); n++){
      if(abs(M(m,n)) != 0){
        cout << "|";
        for(int q = 0; q < qbits; q++){
          cout << (qbit_state(m,qbits,q)?"d":"u");
        }
        cout << "><";
        for(int q = 0; q < qbits; q++){
          cout << (qbit_state(n,qbits,q)?"d":"u");
        }
        cout << "| " << M(m,n) << endl;
      }
    }
  }
}

//--------------------------------------------------------------------------------------------
// NV system objects
//--------------------------------------------------------------------------------------------

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

  spin(Vector3d pos, double g){
    this->pos = pos;
    this->g = g;
  };

  bool operator ==(const spin &s){ return ((pos == s.pos) && (g == s.g)); }
  bool operator !=(const spin &s) { return !(*this == s); }

};

// initialize nitrogen and vacancy centers
const spin n(Vector3d::Zero(),0.);
const spin e(ao,ge);

// unit vectors along lattice directions
const Vector3d zhat = (Vector3d() << 1,1,1).finished()/sqrt(3); // direction from N to V site
const Vector3d xhat = (Vector3d() << 2,-1,-1).finished()/sqrt(6);
const Vector3d yhat = (Vector3d() << 0,1,-1).finished()/sqrt(2);

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// coupling strength between two nuclei; assumes strong magnetic field in zhat
inline double coupling_strength(const spin s1, const spin s2){
  Vector3d r = s2.pos - s1.pos;
  return s1.g*s2.g/(8*pi*pow(r.norm()*a0,3)) * (1 - 3*pow(normed(r).dot(zhat),2));
}

// hyperfine field experienced by nuclear spin s
Vector3d A(const spin s){
  Vector3d r = s.pos - e.pos;
  return e.g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*normed(r).dot(zhat)*normed(r));
}


// group spins into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<spin>> get_clusters(vector<spin> spins, double min_coupling_strength){

  vector<vector<spin>> clusters; // vector of all spin clusters
  vector<int> clustered; // indices of spins we have already clustered

  // loop over indices of all spins
  for(int i = 0; i < spins.size(); i++){
    if(!in_vector(i,clustered)){

      // initialize current cluster
      vector<int> cluster_indices;
      vector<spin> cluster;

      cluster_indices.push_back(i);
      cluster.push_back(spins.at(i));
      clustered.push_back(i);

      // loop over all indices of cluster
      for(vector<int>::size_type ci = 0; ci < cluster.size(); ci++) {

        // loop over all spins indices greater than cluster_indices(ci)
        for(int k = cluster_indices.at(ci)+1; k < spins.size(); k++){

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
int largest_cluster_size(vector<vector<spin>> clusters){
  int largest_size = 0;
  for(int i = 0; i < clusters.size(); i++){
    if(clusters.at(i).size() > largest_size) largest_size = clusters.at(i).size();
  }
  return largest_size;
}

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(vector<spin> spins, int cluster_size_target,
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

// computer NV coherence at scanning frequency w_DD, harmonic k_DD, and fourier component f_DD
double compute_coherence(vector<vector<spin>> clusters,
                         double w_DD, harmonic k_DD, double f_DD, double Bz, int ms){

  double t_DD = 2*pi/w_DD;

  vector<double> ts = pulse_times(k_DD,f_DD);
  double t1 = ts.at(0)*t_DD;
  double t2 = ts.at(1)*t_DD;
  double t3 = t_DD/4;

  MatrixXcd rho_NV_0 = (up+dn)*(up+dn).adjoint()/2;

  double coherence = 1;
  for(int c = 0; c < clusters.size(); c++){
    vector<spin> cluster = clusters.at(c);
    int spins = cluster.size()+1; // total number of spins
    int D = pow(2,spins); // dimensionality of Hilbert space

    // density matrix of initial state
    MatrixXcd rho_0 = act(rho_NV_0, {0}, spins);
    rho_0 /= real(trace(rho_0));

    // construct Hamiltonians
    MatrixXcd H_nn = MatrixXcd::Zero(D,D); // internuclear coupling Hamiltonian
    MatrixXcd H_nZ_eff = MatrixXcd::Zero(D,D); // nuclear Zeeman Hamiltonian
    MatrixXcd H_int = MatrixXcd::Zero(D,D); // interaction Hamiltonian with F(t) = 1

    for(int s = 0; s < cluster.size(); s++){
      for(int r = 0; r < s; r++){
        MatrixXcd H_nn_sr = coupling_strength(cluster.at(s),cluster.at(r))
          * (3*tp(dot(s_vec,zhat),dot(s_vec,zhat)) - dot(s_vec,s_vec));
        H_nn += act(H_nn_sr, {s+1,r+1}, spins);
      }
      const MatrixXcd AI = dot(A(cluster.at(s)), s_vec);
      H_nZ_eff += act( -cluster.at(s).g*dot(Bz*zhat, s_vec) + ms/2*AI, {s+1}, spins);
      H_int += ms/2 * act(tp(sz,AI), {0,s+1}, spins);
    }

    // compute propagator for AXY sequence
    MatrixXcd U1 = exp(-j*(H_nn + H_nZ_eff + H_int) * (t1));
    MatrixXcd U2 = exp(-j*(H_nn + H_nZ_eff - H_int) * (t2-t1));
    MatrixXcd U3 = exp(-j*(H_nn + H_nZ_eff + H_int) * (t3-t2));
    MatrixXcd U4 = exp(-j*(H_nn + H_nZ_eff - H_int) * (t3-t2));
    MatrixXcd U5 = exp(-j*(H_nn + H_nZ_eff + H_int) * (t2-t1));
    MatrixXcd U6 = exp(-j*(H_nn + H_nZ_eff - H_int) * (t1));

    MatrixXcd U = U1*U2*U3*U4*U5*U6 * U6*U5*U4*U3*U2*U1;

    // compute density matrix after a single AXY sequence
    MatrixXcd rho = U*rho_0*U.adjoint();

    // update coherence
    coherence *= pow(2,cluster.size())*real(trace(rho*rho_0));
    cout << pow(2,cluster.size())*real(trace(rho_0*rho_0)) << "     "
         << pow(2,cluster.size())*real(trace(rho*rho)) << "     "
         << pow(2,cluster.size())*real(trace(rho*rho_0)) << endl;
  }
  return coherence;
}
