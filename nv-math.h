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

//--------------------------------------------------------------------------------------------
// Spin vectors and structs
//--------------------------------------------------------------------------------------------

// spin vector for a spin-1/2 particle
const mvec s_vec = mvec(sx,xhat) + mvec(sy,yhat) + mvec(sz,zhat);

// perform spin-1/2 rotation about arbitrary axis
Matrix2cd rotate(const Vector3d axis, const double phi);

// rotate into one axis from another
inline Matrix2cd rotate(const Vector3d axis_end, const Vector3d axis_start){
  return rotate(hat(axis_start.cross(axis_end)), acos(dot(hat(axis_start),hat(axis_end))));
}

// rotate into one basis from another
Matrix2cd rotate(const vector<Vector3d> basis_end, const vector<Vector3d> basis_start);

// struct for spins
struct spin{
  const Vector3d pos; // position
  const double g; // gyromagnetic ratio
  const mvec S; // spin vector

  spin(const Vector3d pos, const double g, const mvec S);

  bool operator==(const spin& s) const {
    return ((pos == s.pos) && (g == s.g) && (S == s.S));
  }
  bool operator!=(const spin& s) const { return !(*this == s); }
};


// harmonic for AXY sequence
enum axy_harmonic { first = 1, third = 3 };

// struct containing system and simulation info
struct nv_system{
  const spin n;
  const spin e;
  const int ms;
  const double static_Bz;
  const double scale_factor;
  const uint integration_factor;

  vector<spin> nuclei;
  double cluster_coupling;
  vector<vector<uint>> clusters;

  nv_system(const int ms, const double static_Bz,
            const double scale_factor, const uint integration_factor);
};

//--------------------------------------------------------------------------------------------
// Spin clustering methods
//--------------------------------------------------------------------------------------------

// determine whether two spins are a larmor pair
bool is_larmor_pair(const nv_system& nv, const uint idx1, const uint idx2);

// coupling strength between two spins; assumes strong magnetic field in zhat
double coupling_strength(const spin& s1, const spin& s2);

// group nuclei into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<uint>> cluster_nuclei(const vector<spin>& nuclei,
                                    const double min_coupling_strength);

// group together clusters close nuclei have similar larmor frequencies
vector<vector<uint>> group_clusters(const nv_system& nv);

// get size of largest spin cluster
uint largest_cluster_size(const vector<vector<uint>>& clusters);

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(const vector<spin>& nuclei, const double initial_cluster_coupling,
                            const uint cluster_size_target, const double dcc_cutoff);

uint get_cluster_containing_index(const nv_system& nv, const uint index);

uint get_index_in_cluster(const uint index, const vector<uint> cluster);

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// hyperfine field experienced by target nucleus
inline Vector3d A(const nv_system& nv, const spin& s){
  const Vector3d r = s.pos - nv.e.pos;
  return nv.e.g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*dot(hat(r),zhat)*hat(r));
};
inline Vector3d A(const nv_system& nv, const uint index){
  return A(nv,nv.nuclei.at(index));
};

// component of hyperfine field perpendicular to the larmor axis
Vector3d A_perp(const nv_system&nv, const spin& s);
inline Vector3d A_perp(const nv_system&nv, const uint index){
  return A_perp(nv,nv.nuclei.at(index));
};

// effective larmor frequency of target nucleus
inline Vector3d effective_larmor(const nv_system& nv, const spin& s){
  return s.g*nv.static_Bz*zhat - nv.ms/2.*A(nv,s);
};
inline Vector3d effective_larmor(const nv_system& nv, const uint index){
  return effective_larmor(nv,nv.nuclei.at(index));
};

// minimum difference in larmor frequencies between target nucleus and other nuclei
//   i.e. min{ |w_s - w_{index}| for all s != index }
double larmor_resolution(const nv_system& nv, const uint index);

// return maximum allowable value of pi*f_k for the AXY sequence
inline double axy_f_max(const axy_harmonic k){ return (k == first ? 8*cos(pi/9)-4 : 4)/pi; };

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const double f, const axy_harmonic k);

// advance pulses by a given (normalized) time
vector<double> advanced_pulse_times(const vector<double> pulse_times, const double advance);

// evaluate F(x) (i.e. sign in front of sigma_z^{NV}) for given AXY pulses
int F_AXY(const double x, const vector<double> pulses);

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2);

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const nv_system& nv, const uint cluster_index);

// NV zero-field splitting plus Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const nv_system& nv, const Vector3d& B){
  return NV_ZFS*dot(nv.e.S,zhat)*dot(nv.e.S,zhat) - nv.e.g*dot(B,nv.e.S);
}

// NV rotation about zhat at a given time
inline MatrixXcd U_NV_GS(const nv_system& nv, const double time, const uint spins = 1){
  return act(exp(-j*time*H_NV_GS(nv,nv.static_Bz*zhat)),{0},spins);
}

// nuclear Zeeman Hamiltonian
MatrixXcd H_nZ(const nv_system& nv, const uint cluster_index, const Vector3d& B);

// Zeeman Hamiltonian for NV center with cluster
inline MatrixXcd H_Z(const nv_system& nv, const uint cluster_index, const Vector3d& B);

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const axy_harmonic k_DD, const double scan_time);

//--------------------------------------------------------------------------------------------
// Control fields and simulation
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

  Vector3d B(const double t) const {
    Vector3d net_B = Vector3d::Zero();
    for(uint c = 0; c < num(); c++){
      net_B += Bs.at(c) * cos(freqs.at(c)*t + phases.at(c));
    }
    return net_B;
  }
};

// return control field for decoupling a single nucleus from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd);

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const double f_DD, const axy_harmonic k_DD,
                              const double simulation_time, const double advance = 0);

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const double f_DD, const axy_harmonic k_DD,
                              const double simulation_time, const control_fields& controls,
                              const double advance = 0);
