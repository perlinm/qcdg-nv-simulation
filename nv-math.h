#pragma once
#define EIGEN_USE_MKL_ALL

#include <eigen3/Eigen/Dense> // linear algebra library

#include "constants.h"
#include "qp-math.h"

using namespace std;
using namespace Eigen;

// ---------------------------------------------------------------------------------------
// Diamond lattice parameters
// ---------------------------------------------------------------------------------------

// diamond lattice vectors in units of the lattice parameter at 300 K
const Vector3d ao = (Vector3d() << 1,1,1).finished()/2;
const Vector3d a1 = (Vector3d() << 0,1,1).finished();
const Vector3d a2 = (Vector3d() << 1,0,1).finished();
const Vector3d a3 = (Vector3d() << 1,1,0).finished();

// unit vectors along bonding axes
const Vector3d zhat = hat(ao); // direction from V to N
const Vector3d xhat = hat(a1-a2);
const Vector3d yhat = zhat.cross(xhat);

// locations of nitrogen nucleus and electron
const Vector3d n_pos = ao;
const Vector3d e_pos = Vector3d::Zero();

// print vec in the {xhat,yhat,zhat} basis
inline Vector3d in_crystal_basis(const Vector3d& vec) {
  return (Vector3d() << dot(vec,xhat), dot(vec,yhat), dot(vec,zhat)).finished();
}

// return "integerized" components of a position vector on z axis or in x-y plane
// these vectors allow for comparing positions without worrying about numerical error
inline Vector3i z_int_pos(const Vector3d& pos) {
  Vector3d p = 6*dot(pos,zhat)*zhat;
  return (Vector3i() << round(p(0)), round(p(1)), round(p(2))).finished();
}
inline Vector3i xy_int_pos(const Vector3d& pos) {
  const Vector3d p = 3*(pos-dot(pos,zhat)*zhat);
  return (Vector3i() << round(p(0)), round(p(1)), round(p(2))).finished();
}

// ---------------------------------------------------------------------------------------
// Spin vectors
// ---------------------------------------------------------------------------------------

// pauli spin vector and spin-1/2 operator
const mvec s_vec = mvec(sx,xhat) + mvec(sy,yhat) + mvec(sz,zhat);
const mvec I_vec = s_vec/2;

// perform spin-1/2 rotation about arbitrary axis
inline Matrix2cd rotate(const Vector3d& axis) {
  if (axis.squaredNorm() > 0) return exp(-j*dot(I_vec,axis));
  else return I2;
}
inline Matrix2cd rotate(const double angle, const Vector3d& axis) {
  return rotate(angle*hat(axis));
}

// rotate into one axis from another
inline Matrix2cd rotate(const Vector3d& axis_end, const Vector3d& axis_start) {
  return rotate(acos(dot(hat(axis_start),hat(axis_end))),
                hat(axis_start.cross(axis_end)));
}

// rotate into one basis from another
Matrix2cd rotate(const vector<Vector3d>& basis_end, const vector<Vector3d>& basis_start);

// ---------------------------------------------------------------------------------------
// NV system struct
// ---------------------------------------------------------------------------------------

// harmonic for AXY sequence
enum axy_harmonic { first = 1, third = 3 };

// struct containing system and simulation info
struct nv_system {
  const vector<Vector3d> nuclei;
  const vector<vector<uint>> clusters;
  const int ms;
  const double static_gBz;
  const axy_harmonic k_DD;
  const double scale_factor;
  const double integration_factor;
  const bool no_nn;

  nv_system(const vector<Vector3d>& nuclei, const vector<vector<uint>>& clusters,
            const int ms, const double static_gBz, const axy_harmonic k_DD,
            const double scale_factor, const double integration_factor, const bool no_nn);

  mvec e_S() const {
    return mvec(sx/sqrt(2),xhat) + mvec(ms*sy/sqrt(2),yhat) + mvec(ms*(sz+I2)/2.,zhat);
  };
};

// ---------------------------------------------------------------------------------------
// Spin placement and clustering methods
// ---------------------------------------------------------------------------------------

// determine whether two nuclei of the same species are a larmor pair
bool is_larmor_pair(const vector<Vector3d>& nuclei, const uint idx1, const uint idx2);
inline bool is_larmor_pair(const nv_system& nv, const uint idx1, const uint idx2) {
  return is_larmor_pair(nv.nuclei,idx1,idx2);
}

// check whether a nucleus is addressable
bool can_address(const vector<Vector3d>& nuclei, const uint target);
inline bool can_address(const nv_system& nv, const uint target) {
  return can_address(nv.nuclei, target);
}

// coupling strength between two C-13 nuclei
inline double coupling_strength(const Vector3d& p1, const Vector3d& p2) {
  return g_C13*g_C13/(4*pi*pow((p2-p1).norm()*a0/2,3));
}
inline double coupling_strength(const vector<Vector3d>& nuclei,
                                const uint idx1, const uint idx2) {
  return coupling_strength(nuclei.at(idx1), nuclei.at(idx2));
}

// coupling strength between two C-13 nuclei with a strong static magnetic field
inline double strong_field_coupling(const Vector3d& p1, const Vector3d& p2) {
  const Vector3d r = p2 - p1;
  return 0.5 * coupling_strength(p1,p2) * abs(1-3*dot(hat(r),zhat));
}
inline double strong_field_coupling(const vector<Vector3d>& nuclei,
                                    const uint idx1, const uint idx2) {
  return strong_field_coupling(nuclei.at(idx1), nuclei.at(idx2));
}
// group nuclei into clusters with intercoupling strengths >= min_coupling_strength
vector<vector<uint>> cluster_with_coupling(const vector<Vector3d>& nuclei,
                                           const double min_coupling_strength,
                                           const bool cluster_by_larmor_frequency = true);

// find largest intercluster coupling strength
double get_cluster_coupling(const vector<Vector3d>& nuclei,
                            const vector<vector<uint>>& clusters);

// get size of largest spin cluster
uint largest_cluster_size(const vector<vector<uint>>& clusters);

// largest internuclear coupling
double largest_coupling(const vector<Vector3d>& nuclei);

// minimum allowable cluster size limit
inline double smallest_possible_cluster_size(const vector<Vector3d>& nuclei) {
  return largest_cluster_size(cluster_with_coupling(nuclei, DBL_MAX));
}

// cluster nuclei and set cluster_coupling for a given maximum cluster size
vector<vector<uint>> cluster_nuclei(const vector<Vector3d>& nuclei,
                                    const uint max_cluster_size,
                                    const bool cluster_by_larmor_frequency,
                                    const double initial_cluster_coupling = 100,
                                    const double cc_resolution = 1e-5);

uint get_cluster_containing_target(const nv_system& nv, const uint index);

uint get_index_in_cluster(const uint index, const vector<uint> cluster);
inline uint get_index_in_cluster(const nv_system& nv, const uint index) {
  const uint cluster = get_cluster_containing_target(nv,index);
  return get_index_in_cluster(index, nv.clusters.at(cluster));
}
inline uint get_index_in_subsystem(const nv_system& nv, const uint index) {
  return get_index_in_cluster(nv, index) + 1;
}

// ---------------------------------------------------------------------------------------
// Hyperfine fields and effective larmor vectors
// ---------------------------------------------------------------------------------------

// hyperfine field experienced by C-13 nucleus located at pos
inline Vector3d hyperfine(const Vector3d& pos) {
  const Vector3d r = pos - e_pos;
  return g_e*g_C13/(4*pi*pow(r.norm()*a0/2,3)) * (zhat - 3*dot(hat(r),zhat)*hat(r));
}
inline Vector3d hyperfine(const nv_system& nv, const uint index) {
  return hyperfine(nv.nuclei.at(index));
}

// effective larmor frequency of target C-13 nucleus located at pos
inline Vector3d effective_larmor(const double static_Bz, const int ms,
                                 const Vector3d& pos) {
  return g_C13*static_Bz*zhat - ms/2.*hyperfine(pos);
}
inline Vector3d effective_larmor(const nv_system& nv, const Vector3d& pos) {
  return nv.static_gBz*zhat - nv.ms/2.*hyperfine(pos);
}
inline Vector3d effective_larmor(const nv_system& nv, const uint index) {
  return effective_larmor(nv,nv.nuclei.at(index));
}

// component of hyperfine field parallel to the larmor axis
inline Vector3d hyperfine_parallel(const double static_Bz, const int ms,
                                   const Vector3d& pos) {
  const Vector3d w_eff_hat = hat(effective_larmor(static_Bz,ms,pos));
  const Vector3d A = hyperfine(pos);
  return dot(A,w_eff_hat)*w_eff_hat;
}

// component of hyperfine field perpendicular to the larmor axis
inline Vector3d hyperfine_perp(const double static_Bz, const int ms,
                               const Vector3d& pos) {
  const Vector3d w_eff_hat = hat(effective_larmor(static_Bz,ms,pos));
  const Vector3d A = hyperfine(pos);
  return A - dot(A,w_eff_hat)*w_eff_hat;
}
inline Vector3d hyperfine_perp(const nv_system&nv, const Vector3d& pos) {
  const Vector3d w_eff_hat = hat(effective_larmor(nv,pos));
  const Vector3d A = hyperfine(pos);
  return A - dot(A,w_eff_hat)*w_eff_hat;
}
inline Vector3d hyperfine_perp(const nv_system&nv, const uint index) {
  return hyperfine_perp(nv,nv.nuclei.at(index));
}

// magnitude of hyperfine field parallel to the NV axis
inline double hyperfine_z(const Vector3d& pos) {
  return dot(hyperfine(pos), zhat);
}

// magnitude of hyperfine field perpendicular to the NV axis
inline double hyperfine_xy(const Vector3d& pos) {
  const Vector3d A = hyperfine(pos);
  return (A - dot(A,zhat)*zhat).norm();
}

// minimum difference in larmor frequencies between target nucleus and other nuclei
//   i.e. min{ |w_s - w_{index}| for all s != index }
double larmor_resolution(const nv_system& nv, const uint index);

// ---------------------------------------------------------------------------------------
// Hamiltonians
// ---------------------------------------------------------------------------------------

// Hamiltoninan coupling two spins
MatrixXcd H_ss(const Vector3d& p1, const double g1, const mvec& S1,
               const Vector3d& p2, const double g2, const mvec& S2);

// Hamiltonian coupling NV electron spin to C-13 nuclear spins in a cluster
MatrixXcd H_en(const nv_system& nv, const uint cluster_index);

// Hamiltonian coupling C-13 nuclear spins in a cluster; acts only on cluster
MatrixXcd H_nn(const nv_system& nv, const uint cluster_index);

// spin-spin coupling Hamiltonian for the entire system
inline MatrixXcd H_int(const nv_system& nv, const uint cluster_index) {
  return H_en(nv, cluster_index) + tp(I2,H_nn(nv, cluster_index));
}

// nuclear Zeeman Hamiltonian; acts only on cluster
MatrixXcd H_nZ(const nv_system& nv, const uint cluster_index, const Vector3d& gB);

// NV Zeeman Hamiltonian
inline MatrixXcd H_NV_Z(const nv_system& nv, const Vector3d& gB) {
  return -g_e/g_C13*dot(gB,nv.e_S());
}

// NV zero-field splitting + static Zeeman Hamiltonian
inline MatrixXcd H_NV_GS(const nv_system& nv) {
  return NV_ZFS*dot(nv.e_S(),zhat)*dot(nv.e_S(),zhat) + H_NV_Z(nv,nv.static_gBz*zhat);
}

// net NV Hamiltonian
inline MatrixXcd H_NV(const nv_system& nv, const Vector3d& gB) {
  return H_NV_GS(nv) + H_NV_Z(nv,gB);
}

// total internal system Hamiltonian
inline MatrixXcd H_sys(const nv_system& nv, const uint cluster) {
  return (H_int(nv,cluster) +
          tp(I2, H_nZ(nv,cluster, nv.static_gBz*zhat)) +
          act(H_NV_GS(nv), {0}, nv.clusters.at(cluster).size()+1));
}

// total control Hamiltonian for the entire system
MatrixXcd H_ctl(const nv_system& nv, const uint cluster, const Vector3d& gB_ctl);

// ---------------------------------------------------------------------------------------
// Control field struct
// ---------------------------------------------------------------------------------------

struct control_fields {
  vector<Vector3d> gBs;
  vector<double> freqs;
  vector<double> phases;

  control_fields(){};
  control_fields(const Vector3d& gB, const double freq = 0, const double phase = 0) {
    this->add(gB,freq,phase);
  }
  control_fields(const vector<Vector3d>& gBs, const vector<double>& freqs,
                 const vector<double>& phases) {
    assert((gBs.size() == freqs.size()) && (gBs.size() == phases.size()));
    this->gBs = gBs;
    this->freqs = freqs;
    this->phases = phases;
  }

  void add(const Vector3d& gB, const double freq = 0, const double phase = 0) {
    gBs.push_back(gB);
    freqs.push_back(freq);
    phases.push_back(phase);
  }
  void add(const control_fields& controls) {
    for (uint i = 0; i < controls.num(); i++) {
      gBs.push_back(controls.gBs.at(i));
      freqs.push_back(controls.freqs.at(i));
      phases.push_back(controls.phases.at(i));
    }
  }
  void remove(const uint i) {
    assert(i < gBs.size());
    gBs.erase(gBs.begin() + i);
    freqs.erase(freqs.begin() + i);
    phases.erase(phases.begin() + i);
  }

  uint num() const { return gBs.size(); }

  Vector3d gB(const double t = 0) const {
    Vector3d net_gB = Vector3d::Zero();
    for (uint c = 0; c < num(); c++) {
      net_gB += gBs.at(c) * cos(freqs.at(c)*t + phases.at(c));
    }
    return net_gB;
  }

  bool all_fields_static() const{
    for (uint c = 0; c < num(); c++) {
      if (freqs.at(c) > 0) return false;
    }
    return true;
  }
};

// ---------------------------------------------------------------------------------------
// Protocol struct
// ---------------------------------------------------------------------------------------

struct protocol {
  MatrixXcd U;
  double time;
  uint pulses;

  protocol(){};
  protocol(const MatrixXcd& U, const double time = 0, const uint pulses = 0) {
    this->U = U;
    this->time = time;
    this->pulses = pulses;
  }

  bool operator==(const protocol& p) const {
    return (time == p.time) && (U == p.U) && (pulses == p.pulses);
  }
  bool operator!=(const protocol& p) const { return !(*this == p); }

  protocol operator*(const protocol& p) const {
    return protocol(U * p.U, time + p.time, pulses + p.pulses);
  }

  protocol adjoint() const {
    return protocol(U.adjoint(), time, pulses);
  }
  protocol pow(const uint n) const {
    return protocol(U.pow(n), time*n, pulses*n);
  }
  protocol pow(const long unsigned int n) const {
    return protocol(U.pow(n), time*n, pulses*n);
  }

  static protocol Identity(const uint D) {
    return protocol(MatrixXcd::Identity(D,D));
  }
};

// methods with protocol objects
inline protocol pow(const protocol& U, const uint n) { return U.pow(n); }
inline protocol pow(const protocol& U, const long unsigned int n) { return U.pow(n); }

inline double protocol_fidelity(const vector<protocol>& P,
                                const vector<uint> system = {}) {
  assert(P.size() == 2);
  return gate_fidelity(P.at(false).U, P.at(true).U, system);
}

// ---------------------------------------------------------------------------------------
// AXY protocol methods
// ---------------------------------------------------------------------------------------

// return maximum allowable value of f_k for the AXY sequence
inline double axy_f_max(const axy_harmonic k) {
  if (k == first) return (8*cos(pi/9) - 4) / pi;
  else return 4/pi; // if k == third
}

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const double f, const axy_harmonic k);

// advance pulses by a given (normalized) time
vector<double> advanced_pulse_times(const vector<double> pulse_times, double advance);

// evaluate F(x) (i.e. sign in front of sigma_z^{NV}) for given AXY pulses
int F_AXY(double x, const vector<double> pulses);

// ---------------------------------------------------------------------------------------
// Simulation methods
// ---------------------------------------------------------------------------------------

// return control field for decoupling a single nucleus from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd);

// perform given rotation on the full NV electron spin
protocol rotate_full_NV(const nv_system& nv, const Vector3d& rotation, const uint spins);

// perform rotation of NV electron spin necessary to generate U_NV
protocol act_NV(const nv_system& nv, const Matrix2cd& U_NV, const uint spins);

// perform given rotation on the reduced 2-level NV electron spin
inline protocol rotate_NV(const nv_system& nv, const double angle, const Vector3d& axis,
                          const uint spins) {
  return act_NV(nv, rotate(angle, axis), spins);
};

// simulate propagator with static control fields
protocol simulate_AXY(const nv_system& nv, const uint cluster,
                      const double w_DD, const double f_DD, const axy_harmonic k_DD,
                      const double simulation_time, const double advance_time = 0,
                      const Vector3d gB_ctl = Vector3d::Zero());

// simulate propagator with dynamic control fields
protocol simulate_AXY(const nv_system& nv, const uint cluster,
                      const double w_DD, const double f_DD, const axy_harmonic k_DD,
                      const control_fields& controls, const double simulation_time,
                      const double advance_time = 0, const double phi_DD = 0);


// perform NV coherence measurement with a static magnetic field on a single cluster
double cluster_coherence(const nv_system& nv, const uint cluster, const double w_scan,
                         const double f_DD, const double scan_time,
                         const Vector3d& gB_ctl = Vector3d::Zero());

// perform NV coherence measurement with a static magnetic field on the entire system
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time,
                             const Vector3d& gB_ctl = Vector3d::Zero());

// perform NV coherence measurement with control fields on the entire system
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time, const control_fields& controls,
                             const double phi_DD = 0,
                             const double precision_factor = 0.05);
