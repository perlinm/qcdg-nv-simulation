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
// Spin state vectors
//--------------------------------------------------------------------------------------------

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

// spin propagators; Ua corresponds to a Hamiltonian # H = h s_a
inline MatrixXcd Ux(const double ht){ return cos(ht)*I2 - j*sin(ht)*sx; }
inline MatrixXcd Uy(const double ht){ return cos(ht)*I2 - j*sin(ht)*sy; }
inline MatrixXcd Uz(const double ht){ return cos(ht)*I2 - j*sin(ht)*sz; }

// rotation operators
inline MatrixXcd Rx(const double phi){ return Ux(phi/2); }
inline MatrixXcd Ry(const double phi){ return Uy(phi/2); }
inline MatrixXcd Rz(const double phi){ return Uz(phi/2); }

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

  nv_system(const int ms, const double static_Bz, const double scale_factor,
            const uint integration_factor);
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

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const uint k, const double f);

// advance pulses by a given phase
vector<double> advanced_pulse_times(const vector<double> pulse_times, const double advance);

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const spin& s1, const spin& s2);

// Hamiltoninan coupling two spins; assumes large static Bz
inline MatrixXcd H_ss_large_static_Bz(const spin& s1, const spin& s2){
  return coupling_strength(s1,s2) * (3*tp(dot(s1.S,zhat), dot(s2.S,zhat)) - dot(s1.S,s2.S));
}

// spin-spin coupling Hamiltonian for NV center with cluster
MatrixXcd H_int(const nv_system& nv, const uint cluster_index);

// spin-spin coupling Hamiltonian for NV center with cluster; assumes large static Bz
MatrixXcd H_int_large_static_Bz(const nv_system& nv, const uint cluster_index);

// NV zero-field splitting plus Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const nv_system& nv, const Vector3d& B);

// nuclear Zeeman Hamiltonian
MatrixXcd H_nZ(const nv_system& nv, const uint cluster_index, const Vector3d& B);

// Zeeman Hamiltonian for NV center with cluster
inline MatrixXcd H_Z(const nv_system& nv, const uint cluster_index, const Vector3d& B);

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const nv_system& nv, const double w_scan, const uint k_DD,
                             const double f_DD, const double scan_time);

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

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const uint k_DD, const double f_DD,
                              const double simulation_time, const double advance = 0);

MatrixXcd simulate_propagator(const nv_system& nv, const uint cluster,
                              const double w_DD, const uint k_DD, const double f_DD,
                              const double simulation_time, const control_fields& controls,
                              const double advance = 0);

//--------------------------------------------------------------------------------------------
// Nuclear targeting methods
//--------------------------------------------------------------------------------------------

// return "natural" basis of a nucleus
vector<Vector3d> natural_basis(const nv_system& nv, const uint index);

// return axis with given azimuth and polar angles in a given basis
inline Vector3d axis(const double azimuth, const double polar = pi/2,
                     const vector<Vector3d> basis = {xhat, yhat, zhat}){
  return cos(polar) * basis.at(2) + sin(polar) * ( cos(azimuth) * basis.at(0) +
                                                   sin(azimuth) * basis.at(1) );
}

// return axis with given polar and azimuthal angles in the "natural" basis
inline Vector3d natural_axis(const nv_system& nv, const uint index,
                             const double azimuth, const double polar = pi/2){
  return axis(azimuth, polar, natural_basis(nv,index));
}

// return control field for decoupling spin s from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd = 0,
                                        const double theta_rfd = pi/2);

// return AXY sequence pulses with given offset
vector<double> offset_pulses(vector<double> xs, const double x_offset);

// propagator generated by t_0 H = phi * sigma_n^{NV}
inline Matrix2cd U_NV(const Vector3d axis, const double phi){
  return exp(-j*phi*dot(s_vec,hat(axis)));
}

// propagator U = exp(-i * rotation_angle * sigma_{axis}^{index})
MatrixXcd U_ctl(const nv_system& nv, const uint index, const double target_axis_azimuth,
                const double rotation_angle, const bool exact = false,
                const bool adjust_AXY = true);

// propagator U = exp(-i * rotation_angle * sigma_{n_1}^{NV}*sigma_{n_2}^{index})
MatrixXcd U_int(const nv_system& nv, const uint index, const uint k_DD,
                const double nv_axis_polar, const double nv_axis_azimuth,
                const double target_axis_azimuth, const double rotation_angle,
                const bool exact = false);

// iSWAP operation
MatrixXcd iSWAP(const nv_system& nv, const uint index, const uint k_DD,
                const bool exact = false);

// SWAP_NVST operation
MatrixXcd SWAP_NVST(const nv_system& nv, const uint idx1, const uint idx2,
                    const uint k_DD, const bool exact = false);

// polarize an arbitrary state into the pure state psi; warning: not a unitary operation
inline MatrixXcd polarize(const VectorXcd psi){
  return psi*VectorXcd::Ones(psi.size()).adjoint();
}
