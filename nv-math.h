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
const Vector3d zhat = hat(ao); // direction from V to N
const Vector3d xhat = hat(a1-a2);
const Vector3d yhat = zhat.cross(xhat);

// print vec in the {xhat,yhat,zhat} basis
inline Vector3d in_basis(const Vector3d& vec){
  return (Vector3d() << dot(vec,xhat), dot(vec,yhat), dot(vec,zhat)).finished();
}

//--------------------------------------------------------------------------------------------
// Spin vectors and structs
//--------------------------------------------------------------------------------------------

// spin vector for a spin-1/2 particle
const mvec s_vec = mvec(sx,xhat) + mvec(sy,yhat) + mvec(sz,zhat);

// perform spin-1/2 rotation about arbitrary axis
inline Matrix2cd rotate(const Vector3d& axis){
  if(axis.squaredNorm() > 0) return exp(-j*dot(s_vec/2.,axis));
  else return I2;
}
inline Matrix2cd rotate(const double angle, const Vector3d& axis){
  return rotate(angle*hat(axis));
}

// rotate into one axis from another
inline Matrix2cd rotate(const Vector3d& axis_end, const Vector3d& axis_start){
  return rotate(acos(dot(hat(axis_start),hat(axis_end))), hat(axis_start.cross(axis_end)));
}

// rotate into one basis from another
Matrix2cd rotate(const vector<Vector3d>& basis_end, const vector<Vector3d>& basis_start);

// struct for spins
struct spin{
  const Vector3d pos; // position
  const double g; // gyromagnetic ratio
  const mvec S; // spin vector

  spin(const Vector3d& pos, const double g, const mvec& S);

  bool operator==(const spin& s) const {
    return ((pos == s.pos) && (g == s.g) && (S == s.S));
  }
  bool operator!=(const spin& s) const { return !(*this == s); }
};


// harmonic for AXY sequence
enum axy_harmonic { first = 1, third = 3 };

// struct containing system and simulation info
struct nv_system{
  const spin n = spin(ao, 0., s_vec/2);
  const spin e;
  const int ms;
  const double static_Bz;
  const axy_harmonic k_DD;
  const double scale_factor;
  const double integration_factor;
  const bool no_nn;

  vector<spin> nuclei;
  double cluster_coupling;
  vector<vector<uint>> clusters;

  nv_system(const int ms, const double static_Bz, const axy_harmonic k_DD,
            const double scale_factor, const double integration_factor, const bool no_nn);
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

// group together clusters sharing larmor pairs
vector<vector<uint>> group_clusters(const nv_system& nv);

// get size of largest spin cluster
uint largest_cluster_size(const vector<vector<uint>>& clusters);

// find cluster coupling for which the largest cluster is >= cluster_size_target
double find_target_coupling(const vector<spin>& nuclei, const double initial_cluster_coupling,
                            const uint cluster_size_target, const double dcc_cutoff);

uint get_cluster_containing_target(const nv_system& nv, const uint index);

uint get_index_in_cluster(const uint index, const vector<uint> cluster);
inline uint get_index_in_cluster(const nv_system& nv, const uint index){
  return get_index_in_cluster(index, nv.clusters.at(get_cluster_containing_target(nv,index)));
}

//--------------------------------------------------------------------------------------------
// AXY scanning methods
//--------------------------------------------------------------------------------------------

// hyperfine field experienced by target nucleus
inline Vector3d hyperfine(const nv_system& nv, const spin& s){
  const Vector3d r = s.pos - nv.e.pos;
  return nv.e.g*s.g/(4*pi*pow(r.norm()*a0,3)) * (zhat - 3*dot(hat(r),zhat)*hat(r));
}
inline Vector3d hyperfine(const nv_system& nv, const uint index){
  return hyperfine(nv,nv.nuclei.at(index));
}

// component of hyperfine field perpendicular to the larmor axis
Vector3d hyperfine_perp(const nv_system&nv, const spin& s);
inline Vector3d hyperfine_perp(const nv_system&nv, const uint index){
  return hyperfine_perp(nv,nv.nuclei.at(index));
}

// effective larmor frequency of target nucleus
inline Vector3d effective_larmor(const nv_system& nv, const spin& s){
  return s.g*nv.static_Bz*zhat - nv.ms/2.*hyperfine(nv,s);
}
inline Vector3d effective_larmor(const nv_system& nv, const uint index){
  return effective_larmor(nv,nv.nuclei.at(index));
}

// minimum difference in larmor frequencies between target nucleus and other nuclei
//   i.e. min{ |w_s - w_{index}| for all s != index }
double larmor_resolution(const nv_system& nv, const uint index);

// return maximum allowable value of f_k for the AXY sequence
inline double axy_f_max(const axy_harmonic k){
  if(k == first) return (8*cos(pi/9) - 4) / pi;
  else return 4/pi;
}

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const double f, const axy_harmonic k);

// advance pulses by a given (normalized) time
vector<double> advanced_pulse_times(const vector<double> pulse_times, const double advance);

// evaluate F(x) (i.e. sign in front of sigma_z^{NV}) for given AXY pulses
int F_AXY(const double x, const vector<double> pulses);

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2);

// spin coupling Hamiltonian for the entire spin system
MatrixXcd H_int(const nv_system& nv, const uint cluster_index);

// nuclear Zeeman Hamiltonian
MatrixXcd H_nZ(const nv_system& nv, const uint cluster_index, const Vector3d& B);

// NV Zeeman Hamiltonian
inline MatrixXcd H_NV_Z(const nv_system& nv, const Vector3d& B){
  return -nv.e.g*dot(B,nv.e.S);
}

// NV zero-field splitting + static Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const nv_system& nv){
  return NV_ZFS*dot(nv.e.S,zhat)*dot(nv.e.S,zhat) + H_NV_Z(nv,nv.static_Bz*zhat);
}

// net NV Hamiltonian
inline MatrixXcd H_NV(const nv_system& nv, const Vector3d& B){
  return H_NV_GS(nv) + H_NV_Z(nv,B);
}

// total internal system Hamiltonian
inline MatrixXcd H_sys(const nv_system& nv, const uint cluster){
  return (H_int(nv,cluster) + H_nZ(nv,cluster,nv.static_Bz*zhat) +
          act(H_NV_GS(nv), {0}, nv.clusters.at(cluster).size()+1));
}

// total control Hamiltonian for the entire system
inline MatrixXcd H_ctl(const nv_system& nv, const uint cluster, const Vector3d& B_ctl){
  return H_nZ(nv,cluster,B_ctl) + act(H_NV_Z(nv,B_ctl),{0},nv.clusters.at(cluster).size()+1);
}

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time);

//--------------------------------------------------------------------------------------------
// Control fields and simulation
//--------------------------------------------------------------------------------------------

struct control_fields{
  vector<Vector3d> Bs;
  vector<double> freqs;
  vector<double> phases;

  control_fields(){};
  control_fields(const Vector3d& B, const double freq = 0, const double phase = 0){
    this->add(B,freq,phase);
  }
  control_fields(const vector<Vector3d>& Bs, const vector<double>& freqs,
                 const vector<double>& phases){
    assert((Bs.size() == freqs.size()) && (Bs.size() == phases.size()));
    this->Bs = Bs;
    this->freqs = freqs;
    this->phases = phases;
  }

  void add(const Vector3d& B, const double freq, const double phase = 0){
    Bs.push_back(B);
    freqs.push_back(freq);
    phases.push_back(phase);
  }
  void add(const control_fields& controls){
    for(uint i = 0; i < controls.num(); i++){
      Bs.push_back(controls.Bs.at(i));
      freqs.push_back(controls.freqs.at(i));
      phases.push_back(controls.phases.at(i));
    }
  }
  void remove(const uint i){
    assert(i < Bs.size());
    Bs.erase(Bs.begin() + i);
    freqs.erase(freqs.begin() + i);
    phases.erase(phases.begin() + i);
  }

  uint num() const { return Bs.size(); }

  Vector3d B(const double t) const {
    Vector3d net_B = Vector3d::Zero();
    for(uint c = 0; c < num(); c++){
      net_B += Bs.at(c) * cos(freqs.at(c)*t + phases.at(c));
    }
    return net_B;
  }

  bool all_fields_static() const{
    for(uint c = 0; c < num(); c++){
      if(freqs.at(c) > 0) return false;
    }
    return true;
  }
};

// return control field for decoupling a single nucleus from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd);

// perform given rotation on the NV center
MatrixXcd rotate_NV(const nv_system& nv, const Vector3d& rotation, const uint spins);

// compute and perform rotation of NV center necessary to generate U
MatrixXcd act_NV(const nv_system& nv, const Matrix2cd& U, const uint spins);

// simulate propagator with static control fields
MatrixXcd simulate_AXY8(const nv_system& nv, const uint cluster,
                        const double w_DD, const double f_DD, const axy_harmonic k_DD,
                        const double simulation_time, const double advance = 0,
                        const Vector3d B_ctl = Vector3d::Zero());

// simulate propagator with dynamic control fields
MatrixXcd simulate_AXY8(const nv_system& nv, const uint cluster,
                        const double w_DD, const double f_DD, const axy_harmonic k_DD,
                        const control_fields& controls, const double simulation_time,
                        const double advance = 0);

//--------------------------------------------------------------------------------------------
// Protocol object
//--------------------------------------------------------------------------------------------

struct protocol{
  MatrixXcd U;
  double t;

  protocol(){};
  protocol(const MatrixXcd& U, const double t){
    this->U = U;
    this->t = t;
  }

  bool operator==(const protocol& p) const {
    return (t == p.t) && (U == p.U);
  }
  bool operator!=(const protocol& p) const { return !(*this == p); }
  protocol operator*(const protocol& p) const {
    return protocol(U * p.U, t + p.t);
  }

  protocol adjoint() const {
    return protocol(U.adjoint(), t);
  }
};
