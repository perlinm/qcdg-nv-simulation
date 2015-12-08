#pragma once

using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

//--------------------------------------------------------------------------------------------
// Single spin state vectors
//--------------------------------------------------------------------------------------------

// spin up/down state vectors
const VectorXcd up = (Vector2cd() << 1,0).finished();
const VectorXcd dn = (Vector2cd() << 0,1).finished();

// two qbit basis states
const VectorXcd uu = tp(up,up);
const VectorXcd ud = tp(up,dn);
const VectorXcd du = tp(dn,up);
const VectorXcd dd = tp(dn,dn);

// two qbit singlet-triplet states
const VectorXcd S = (ud-du)/sqrt(2);
const VectorXcd T = (ud+du)/sqrt(2);

//--------------------------------------------------------------------------------------------
// Single spin operations
//--------------------------------------------------------------------------------------------

// pauli spin matrices
const MatrixXcd st = up*up.adjoint() + dn*dn.adjoint();
const MatrixXcd sx = up*dn.adjoint() + dn*up.adjoint();
const MatrixXcd sy = j*(-up*dn.adjoint() + dn*up.adjoint());
const MatrixXcd sz = up*up.adjoint() - dn*dn.adjoint();

// spin propagators; Ua corresponds to a Hamiltonian # H = h s_a
inline MatrixXcd Ux(double ht){ return cos(ht)*I2 - j*sin(ht)*sx; }
inline MatrixXcd Uy(double ht){ return cos(ht)*I2 - j*sin(ht)*sy; }
inline MatrixXcd Uz(double ht){ return cos(ht)*I2 - j*sin(ht)*sz; }

// rotation operators
inline MatrixXcd Rx(double phi){ return Ux(phi/2); }
inline MatrixXcd Ry(double phi){ return Uy(phi/2); }
inline MatrixXcd Rz(double phi){ return Uz(phi/2); }

//--------------------------------------------------------------------------------------------
// Diamond lattice parameters
//--------------------------------------------------------------------------------------------

// diamond lattice vectors, scaled to the lattice parameter (i.e. the unit cell side length)
const Vector3d ao = (Vector3d() << 1,1,1).finished()/4;
const Vector3d a1 = (Vector3d() << 0,1,1).finished()/2;
const Vector3d a2 = (Vector3d() << 1,0,1).finished()/2;
const Vector3d a3 = (Vector3d() << 1,1,0).finished()/2;

// vector of lattice sites in a diamond unit cell
const vector<Vector3d> cell_sites { Vector3d::Zero(), a1, a2, a3, ao, ao+a1, ao+a2, ao+a3 };

// unit vectors along bonding axes
const Vector3d zhat = (Vector3d() << 1,1,1).finished()/sqrt(3); // direction from V to N
const Vector3d xhat = (Vector3d() << 2,-1,-1).finished()/sqrt(6);
const Vector3d yhat = (Vector3d() << 0,1,-1).finished()/sqrt(2);

//--------------------------------------------------------------------------------------------
// Spin vectors and structs
//--------------------------------------------------------------------------------------------

// spin vector for a spin-1/2 particle
const mvec s_vec = mvec(sx/2,xhat) + mvec(sy/2,yhat) + mvec(sz/2,zhat);

// struct for spins
struct spin{
  Vector3d pos; // position
  double g; // gyromagnetic ratio
  mvec S; // spin vector

  spin(const Vector3d& pos, double g, const mvec& S){
    assert(S.size() == 3);
    this->pos = pos;
    this->g = g;
    this->S = S;
  };

  bool operator==(const spin& s) const {
    return ((pos == s.pos) && (g == s.g) && (S == s.S));
  }
  bool operator!=(const spin& s) const { return !(*this == s); }

};

// initialize nitrogen and vacancy centers
const spin n(ao, 0., s_vec);
inline spin e(const int ms){
  return spin(Vector3d::Zero(), ge,
              mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat));
}

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// coupling strength between two spins; assumes strong magnetic field in zhat
double coupling_strength(const spin& s1, const spin& s2);

// group nuclei into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<uint>> get_index_clusters(const vector<spin>& nuclei,
                                        const double min_coupling_strength);

// get size of largest spin cluster
uint largest_cluster_size(const vector<vector<uint>>& ind_clusters);

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(const vector<spin>& nuclei, const uint cluster_size_target,
                            const double initial_cluster_coupling, const double dcc_cutoff);

// group nuclei into clusters according to cluster_indices
vector<vector<spin>> group_nuclei(const vector<spin>& nuclei,
                                 const vector<vector<uint>>& cluster_indices);

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// hyperfine field experienced by spin s
inline Vector3d A(const spin& s){
  const Vector3d r = s.pos - e(1).pos;
  return e(1).g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*dot(hat(r),zhat)*hat(r));
};

// effective larmor frequency of spin s
inline Vector3d effective_larmor(const spin& s, const double static_B, const int ms){
  return s.g*static_B*zhat - ms/2.*A(s);
};

// pulse times for harmonic h and fourier component f
vector<double> pulse_times(const uint k, const double f);

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2);

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const spin& e, const vector<spin>& cluster);

// NV zero-field splitting plus Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const spin& e, const vector<spin>& cluster, const Vector3d& B);

// nuclear Zeeman Hamiltonian
MatrixXcd H_nZ(const vector<spin>& cluster, const Vector3d& B);

// Zeeman Hamiltonian for NV center with cluster
inline MatrixXcd H_Z(const spin& e, const vector<spin>& cluster, const Vector3d& B);

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const double scan_time, const vector<vector<spin>>& clusters,
                             const double w_scan, const uint k_DD, const double f_DD,
                             const double static_B, const int ms);

//--------------------------------------------------------------------------------------------
// Control field scanning and targeting
// (assume large static magnetic field along NV axis)
//--------------------------------------------------------------------------------------------

struct control_fields{
  vector<Vector3d> Bs;
  vector<double> freqs;
  vector<double> phases;

  control_fields(){};
  control_fields(const Vector3d& B, const double w, const double phase){
    this->add(B,w,phase);
  };
  control_fields(const vector<Vector3d>& Bs, const vector<double>& freqs,
                 const vector<double>& phases){
    assert((Bs.size() == freqs.size()) && (Bs.size() == phases.size()));
    this->Bs = Bs;
    this->freqs = freqs;
    this->phases = phases;
  };

  void add(const Vector3d& B, const double w, const double phase){
    Bs.push_back(B);
    freqs.push_back(w);
    phases.push_back(phase);
  };
  void add(const control_fields& controls){
    for(uint i = 0; i < controls.num(); i++){
      Bs.push_back(controls.Bs.at(i));
      freqs.push_back(controls.freqs.at(i));
      phases.push_back(controls.phases.at(i));
    }
  };
  void remove(const uint i){
    assert(i < Bs.size());
    Bs.erase(Bs.begin() + i);
    freqs.erase(freqs.begin() + i);
    phases.erase(phases.begin() + i);
  };

  uint num() const { return Bs.size(); }
};

// Hamiltoninan coupling two spins
inline MatrixXcd H_ss_large_static_B(const spin& s1, const spin& s2);

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int_large_static_B(const spin& e, const vector<spin>& cluster);

// perform NV coherence measurement with a static magnetic field and additional control fields
double coherence_measurement(double scan_time, const vector<vector<spin>>& clusters,
                             const double w_scan, const uint k_DD, const double f_DD,
                             const double static_B, const int ms,
                             const control_fields& controls,
                             const uint integration_factor = 100);

//--------------------------------------------------------------------------------------------
// Single nuclear targeting
//--------------------------------------------------------------------------------------------

// return control field for decoupling spin s from other nuclei
control_fields nuclear_decoupling_field(const spin& s, const double static_B, const int ms,
                                        const double phi_rfd = 0,
                                        const double theta_rfd = pi/2,
                                        const double scale = 1000);

// return AXY sequence pulses with given offset
vector<double> offset_pulses(vector<double> xs, const double x_offset);

// compute fidelity of SWAP operation between NV center and target nucleus
double iswap_fidelity(const uint target_nucleus_index, const vector<spin>& nuclei,
                      const double static_B, const int ms,
                      const uint k_DD = 1, const double scale = 1000);
