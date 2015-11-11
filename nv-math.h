#pragma once

using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

//--------------------------------------------------------------------------------------------
// Diamond lattice parameters
//--------------------------------------------------------------------------------------------

// diamond lattice vectors, scaled to the lattice parameter (i.e. the unit cell side length)
const Vector3d ao = (Vector3d() << 1,1,1).finished()/4;
const Vector3d a1 = (Vector3d() << 0,1,1).finished()/2;
const Vector3d a2 = (Vector3d() << 1,0,1).finished()/2;
const Vector3d a3 = (Vector3d() << 1,1,0).finished()/2;

// unit vectors along bonding axes
const Vector3d zhat = (Vector3d() << 1,1,1).finished()/sqrt(3); // direction from V to N
const Vector3d xhat = (Vector3d() << 2,-1,-1).finished()/sqrt(6);
const Vector3d yhat = (Vector3d() << 0,1,-1).finished()/sqrt(2);

// vector of lattice sites in a diamond unit cell
const vector<Vector3d> cell_sites { Vector3d::Zero(), a1, a2, a3, ao, ao+a1, ao+a2, ao+a3 };

//--------------------------------------------------------------------------------------------
// Spin states, matrices, vectors, and structs
//--------------------------------------------------------------------------------------------

// spin up/down state vectors
const VectorXcd up = (Vector2cd() << 1,0).finished();
const VectorXcd dn = (Vector2cd() << 0,1).finished();

// pauli spin matrices
const MatrixXcd st = up*up.adjoint() + dn*dn.adjoint();
const MatrixXcd sx = up*dn.adjoint() + dn*up.adjoint();
const MatrixXcd sy = j*(-up*dn.adjoint() + dn*up.adjoint());
const MatrixXcd sz = up*up.adjoint() - dn*dn.adjoint();

// spin vector for a spin-1/2 particle
const mvec s_vec = mvec(sx/2,xhat) + mvec(sy/2,yhat) + mvec(sz/2,zhat);

// struct for spins
struct spin{
  Vector3d pos; // position
  double g; // gyromagnetic ratio
  mvec S; // spin vector

  spin(const Vector3d pos, double g, const mvec S){
    assert(S.size() == 3);
    this->pos = pos;
    this->g = g;
    this->S = S;
  };

  bool operator==(const spin& s) const {
    return ((pos == s.pos)& & (g == s.g)& & (S == s.S));
  }
  bool operator!=(const spin& s) const { return !(*this == s); }

};

// initialize nitrogen and vacancy centers
const spin n(ao, 0., s_vec);
inline spin e(int ms){
  return spin(Vector3d::Zero(), ge,
              mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat));
}

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// coupling strength between two spins; assumes strong magnetic field in zhat
double coupling_strength(const spin s1, const spin s2);

// group spins into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<spin>> get_clusters(vector<spin> spins, double min_coupling_strength);

// get size of largest spin cluster
uint largest_cluster_size(vector<vector<spin>> clusters);

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(vector<spin> spins, uint cluster_size_target,
                            double initial_cluster_coupling, double dcc_cutoff);

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// hyperfine field experienced by spin s
Vector3d A(const spin s);

// effective larmor frequency of spin s
double effective_larmor(const spin s, const Vector3d B, const Vector3d A, int ms);

// pulse times for harmonic h and fourier component f
vector<double> pulse_times(uint k, double f);

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin s1, const spin s2);

// return spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const spin e, const vector<spin> cluster);

// return Zeeman Hamiltonian for NV center with cluster
MatrixXcd H_Z(const spin e, const vector<spin> cluster, Vector3d B);

// perform NV coherence measurement
double coherence_measurement(int ms, vector<vector<spin>> clusters, double w_scan,
                             uint k_DD, double f_DD, double Bz, double scan_time);
