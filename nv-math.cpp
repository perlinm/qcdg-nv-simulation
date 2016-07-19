#define EIGEN_USE_MKL_ALL

#include <iostream> // for standard output
#include <eigen3/Eigen/Dense> // linear algebra library

#include "constants.h"
#include "qp-math.h"
#include "nv-math.h"

using namespace std;
using namespace Eigen;

// ---------------------------------------------------------------------------------------
// Spin vectors and structs
// ---------------------------------------------------------------------------------------

// rotate into one basis from another
Matrix2cd rotate(const vector<Vector3d>& basis_end, const vector<Vector3d>& basis_start) {
  assert(basis_start.size() == 3);
  assert(basis_end.size() == 3);
  assert(dot(basis_start.at(0).cross(basis_start.at(1)),basis_start.at(2)) > 0);
  assert(dot(basis_end.at(0).cross(basis_end.at(1)),basis_end.at(2)) > 0);

  // rotation matrix taking starting basis vectors to ending basis vectors
  const Matrix3d rotation = (basis_end.at(0)*basis_start.at(0).transpose() +
                             basis_end.at(1)*basis_start.at(1).transpose() +
                             basis_end.at(2)*basis_start.at(2).transpose());
  // rotation angle
  const double angle = acos((rotation.trace()-1.)/2.);

  // determine rotation axis
  const double axis_x = rotation(2,1)-rotation(1,2);
  const double axis_y = rotation(0,2)-rotation(2,0);
  const double axis_z = rotation(1,0)-rotation(0,1);
  const Vector3d axis = hat((Vector3d() << axis_x, axis_y, axis_z).finished());
  if (axis.squaredNorm() > 0) {
    return rotate(angle, axis);
  } else {
    const EigenSolver<Matrix3d> solver(rotation);
    const Vector3cd e_vals = solver.eigenvalues();
    const Matrix3cd e_vecs = solver.eigenvectors();
    uint axis_index = 0;
    for (uint i = 1; i < e_vals.size(); i++) {
      if (abs(e_vals(i)-1.) < abs(e_vals(axis_index)-1.)) axis_index = i;
    }
    return rotate(angle, e_vecs.col(axis_index).real());
  }
}

nv_system::nv_system(const vector<Vector3d>& nuclei, const vector<vector<uint>>& clusters,
                     const int ms, const double static_gBz, const axy_harmonic k_DD,
                     const double scale_factor, const double integration_factor,
                     const bool no_nn) :
  nuclei(nuclei), clusters(clusters), ms(ms), static_gBz(static_gBz), k_DD(k_DD),
  scale_factor(scale_factor), integration_factor(integration_factor), no_nn(no_nn)
{};

// ---------------------------------------------------------------------------------------
// Spin placement and clustering methods
// ---------------------------------------------------------------------------------------

// determine whether two nuclei of the same species are a larmor pair
bool is_larmor_pair(const vector<Vector3d>& nuclei, const uint idx1, const uint idx2) {
  const Vector3d r1 = nuclei.at(idx1) - e_pos;
  const Vector3d r2 = nuclei.at(idx2) - e_pos;
  return
    z_int_pos(r1).squaredNorm() == z_int_pos(r2).squaredNorm() &&
    xy_int_pos(r1).squaredNorm() == xy_int_pos(r2).squaredNorm();
}

// check whether a nucleus is addressable
bool can_address(const vector<Vector3d>& nuclei, const uint target) {
  const Vector3d r = nuclei.at(target) - e_pos;
  const Vector3i r_xy = xy_int_pos(r);
  const Vector3i r_z = z_int_pos(r);

  // exclude targets on the z axis or in the x-y plane
  if (r_xy.squaredNorm() == 0 || r_z.squaredNorm() == 0) return false;

  for (uint i = 0; i < nuclei.size(); i++) {
    if (i == target) continue;

    if (is_larmor_pair(nuclei,target,i)) {
      // exclude larmor pairs with parallel x-y components of the hyperfine field
      const Vector3d s = nuclei.at(i) - e_pos;
      const Vector3i s_xy = xy_int_pos(s);
      if (s_xy == r_xy || s_xy == -r_xy) return false;

      // exclude larmor triplets
      for (uint j = i+1; j < nuclei.size(); j++) {
        if (j == target) continue;
        if (is_larmor_pair(nuclei,i,j)) return false;
      }
    }
  }

  return true;
}

// group nuclei into clusters with intercoupling strengths >= coupling_strength
vector<vector<uint>> cluster_with_coupling(const vector<Vector3d>& nuclei,
                                           const double min_coupling_strength,
                                           const bool cluster_larmor_pairs) {

  vector<vector<uint>> clusters(0); // all clusters
  vector<uint> clustered(0); // indices of nuclei we have already clustered

  for (uint n_seed = 0; n_seed < nuclei.size(); n_seed++) {
    if (in_vector(n_seed,clustered)) continue;

    // initialize current cluster, and add the current nucleus
    vector<uint> cluster;
    cluster.push_back(n_seed);
    clustered.push_back(n_seed);

    // loop over all nuclei in cluster
    for (uint i = 0; i < cluster.size(); i++) {
      const uint n_old = cluster.at(i);
      // loop over all nuclei, checking whether to add them to the cluster
      for (uint n_new = 0; n_new < nuclei.size(); n_new++) {
        if (in_vector(n_new,clustered)) continue;

        // if n_old and n_new are interacting or are a larmor pair,
        //   add n_new to this cluster
        if ((coupling_strength(nuclei,n_old,n_new) > min_coupling_strength) ||
             (cluster_larmor_pairs && is_larmor_pair(nuclei,n_old,n_new))) {
          cluster.push_back(n_new);
          clustered.push_back(n_new);
        }
      }
    }
    clusters.push_back(cluster);
  }
  return clusters;
}

// get size of largest cluster
uint largest_cluster_size(const vector<vector<uint>>& clusters) {
  uint largest_size = 0;
  for (uint i = 0; i < clusters.size(); i++) {
    if (clusters.at(i).size() > largest_size) largest_size = clusters.at(i).size();
  }
  return largest_size;
}

// largest internuclear coupling
double largest_coupling(const vector<Vector3d>& nuclei) {
    double max = 0;
    for (uint i = 0; i < nuclei.size(); i++) {
      for (uint j = i+1; j < nuclei.size(); j++) {
        double c_ij = coupling_strength(nuclei,i,j);
        if (c_ij > max) max = c_ij;
      }
    }
    return max;
}

// cluster nuclei and set cluster_coupling for a given maximum cluster size
vector<vector<uint>> cluster_nuclei(const vector<Vector3d>& nuclei,
                                    const uint max_cluster_size,
                                    const bool cluster_larmor_pairs,
                                    const double initial_cluster_coupling,
                                    const double cc_resolution) {
  assert(cc_resolution > 0);
  uint cluster_size_target = max_cluster_size;
  vector<vector<uint>> clusters;

  // special cases for small clusters
  if (cluster_size_target >= nuclei.size()) {
    clusters.push_back({});
    for (uint n = 0; n < nuclei.size(); n++) {
      clusters.at(0).push_back(n);
    }
    return clusters;
  }
  const uint min_cluster_size_cap = smallest_possible_cluster_size(nuclei);
  if (cluster_size_target < min_cluster_size_cap) {
    for (uint n = 0; n < nuclei.size(); n++) {
      clusters.push_back({n});
    }
    return clusters;
  }

  // find largest coupling for which the largest cluster size is <= cluster_size_target
  double cluster_coupling = initial_cluster_coupling;
  double dcc = cluster_coupling/8;

  bool coupling_too_small = true;
  bool last_coupling_too_small = true;
  bool crossed_correct_coupling = false;

  while (dcc >= cc_resolution) {
    if (coupling_too_small == last_coupling_too_small) {
      if (!crossed_correct_coupling) dcc *= 2;
    } else {
      dcc /= 2;
    }
    last_coupling_too_small = coupling_too_small;

    cluster_coupling += coupling_too_small ? dcc : -dcc;
    clusters = cluster_with_coupling(nuclei, cluster_coupling, cluster_larmor_pairs);

    coupling_too_small = (largest_cluster_size(clusters) > cluster_size_target);
    if (!crossed_correct_coupling && (coupling_too_small != last_coupling_too_small)) {
      crossed_correct_coupling = true;
    }
  }

  while (largest_cluster_size(clusters) > cluster_size_target) {
    cluster_coupling += dcc;
    clusters = cluster_with_coupling(nuclei, cluster_coupling, cluster_larmor_pairs);
  }
  return clusters;
}

// find largest intercluster coupling strength
double get_cluster_coupling(const vector<Vector3d>& nuclei,
                            const vector<vector<uint>>& clusters) {
  double max_coupling = 0;
  for (uint c1 = 0; c1 < clusters.size(); c1++) {
    for (uint n1 = 0; n1 < clusters.at(c1).size(); n1++) {
      for (uint c2 = c1+1; c2 < clusters.size(); c2++) {
        for (uint n2 = 0; n2 < clusters.at(c2).size(); n2++) {
          max_coupling = max(max_coupling,
                             coupling_strength(nuclei,
                                               clusters.at(c1).at(n1),
                                               clusters.at(c2).at(n2)));
        }
      }
    }
  }
  return max_coupling;
}

uint get_cluster_containing_target(const nv_system& nv, const uint index) {
  assert(index < nv.nuclei.size());
  for (uint c = 0; c < nv.clusters.size(); c++) {
    for (uint s = 0; s < nv.clusters.at(c).size(); s++) {
      if (nv.clusters.at(c).at(s) == index) {
        return c;
      }
    }
  }
  return 0;
}

uint get_index_in_cluster(const uint index, const vector<uint> cluster) {
  assert(in_vector(index,cluster));
  for (uint s = 0; s < cluster.size(); s++) {
    if (cluster.at(s) == index) {
      return s;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------------------
// AXY scanning methods
// ---------------------------------------------------------------------------------------

// component of hyperfine field perpendicular to the larmor axis
Vector3d hyperfine_perp(const nv_system&nv, const Vector3d& pos) {
  const Vector3d w_eff_hat = hat(effective_larmor(nv,pos));
  const Vector3d A = hyperfine(pos);
  return A - dot(A,w_eff_hat)*w_eff_hat;
}

// minimum difference in larmor frequencies between target nucleus and other nuclei
//   i.e. min{ |w_s - w_{target}| for all s != index }
double larmor_resolution(const nv_system& nv, const uint target) {
  const double target_larmor = effective_larmor(nv,target).norm();
  double dw_min = target_larmor; // maximum allowable larmor resolution
  for (uint s = 0; s < nv.nuclei.size(); s++) {
    if (is_larmor_pair(nv,s,target)) continue; // find dw_min for distinct frequencies only
    const double dw = abs(target_larmor - effective_larmor(nv,s).norm());
    if (dw < dw_min) dw_min = dw;
  }
  return dw_min;
}

// pulse times for harmonic h and fourier component f
vector<double> axy_pulse_times(const double f, const axy_harmonic k) {
  assert(abs(f) <= axy_f_max(k));
  const double fp = f*pi;

  // compute first two pulse times
  double x1,x2;
  if (k == 1) {
    const double w1 = 4 - fp;
    const double w2 = w1 * (960 - 144*fp - 12*fp*fp + fp*fp*fp);

    x1 = 1/(2*pi) * atan2( (3*fp-12)*w1 + sqrt(3*w2),
                           sqrt(6)*sqrt(w2 - 96*fp*w1 + w1*w1*sqrt(3*w2)));
    x2 = 1/(2*pi) * atan2(-(3*fp-12)*w1 + sqrt(3*w2),
                          sqrt(6)*sqrt(w2 - 96*fp*w1 - w1*w1*sqrt(3*w2)));
  } else { // if k == 3
    const double q1 = 4/(sqrt(5+fp)-1);
    const double q2 = 4/(sqrt(5+fp)+1);

    x1 = 1./4 - 1/(2*pi)*atan(sqrt(q1*q1-1));
    x2 = 1./4 - 1/(2*pi)*atan(sqrt(q2*q2-1));
  }

  // construct vector of all pulse times (normalized to one AXY period)
  vector<double> pulse_times;
  pulse_times.push_back(0);
  pulse_times.push_back(x1);
  pulse_times.push_back(x2);
  pulse_times.push_back(0.25);
  pulse_times.push_back(0.5 - x2);
  pulse_times.push_back(0.5 - x1);
  pulse_times.push_back(0.5 + x1);
  pulse_times.push_back(0.5 + x2);
  pulse_times.push_back(0.75);
  pulse_times.push_back(1. - x2);
  pulse_times.push_back(1. - x1);
  pulse_times.push_back(1);
  return pulse_times;
}

vector<double> advanced_pulse_times(const vector<double> pulse_times, double advance) {
  advance = mod(advance);
  if (advance == 0) return pulse_times;

  // number of pulses
  const uint N = pulse_times.size()-2;

  // advanced pulse_times
  vector<double> advanced_pulse_times;
  advanced_pulse_times.push_back(0);
  for (uint p = 0; p < 2*N; p++) {
    if (p/N + pulse_times.at(p%N+1) - advance >= 0) {
      advanced_pulse_times.push_back(p/N + pulse_times.at(p%N+1) - advance);
    }
    if (advanced_pulse_times.size()-1 == N) break;
  }
  advanced_pulse_times.push_back(1);
  return advanced_pulse_times;
}

// evaluate F(x) (i.e. sign in front of sigma_z^{NV}) for given AXY pulses
int F_AXY(double x, const vector<double> pulses) {
  x = mod(x);
  uint pulse_count = 0;
  for (uint i = 1; i < pulses.size()-1; i++) {
    if (pulses.at(i) < x) pulse_count++;
    else break;
  }
  return (pulse_count % 2 == 0 ? 1 : -1);
}

// Hamiltoninan coupling two C-13 nuclei
MatrixXcd H_ss(const Vector3d& p1, const double g1, const mvec& S1,
               const Vector3d& p2, const double g2, const mvec& S2) {
  const Vector3d r = p2 - p1;
  return g1*g2/(4*pi*pow(r.norm()*a0/2,3))
    * (dot(S1,S2) - 3*tp(dot(S1,hat(r)), dot(S2,hat(r))));
}

// Hamiltonian coupling NV electron to C-13 nuclei in a cluster
MatrixXcd H_en(const nv_system& nv, const uint cluster_index) {
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  const int spins = cluster.size()+1;
  MatrixXcd H = MatrixXcd::Zero(pow(2,spins),pow(2,spins));
  for (uint n1 = 0; n1 < cluster.size(); n1++) {
    const Vector3d n1_pos = nv.nuclei.at(cluster.at(n1));
    H += act(H_ss(e_pos, g_e, nv.e_S(), n1_pos, g_C13, I_vec),
             {0,n1+1}, spins);
  }
  return H;
}

// Hamiltonian coupling C-13 nuclei in a cluster; acts only on cluster
MatrixXcd H_nn(const nv_system& nv, const uint cluster_index){
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  MatrixXcd H = MatrixXcd::Zero(pow(2,cluster.size()),pow(2,cluster.size()));
  for (uint n1 = 0; n1 < cluster.size(); n1++) {
    const Vector3d n1_pos = nv.nuclei.at(cluster.at(n1));
    for (uint n2 = 0; n2 < n1; n2++) {
      const Vector3d n2_pos = nv.nuclei.at(cluster.at(n2));
      H += act(H_ss(n1_pos, g_C13, I_vec, n2_pos, g_C13, I_vec),
               {n1,n2}, cluster.size());
    }
  }
  return H;
}

// spin-spin coupling Hamiltonian for the entire system
MatrixXcd H_int(const nv_system& nv, const uint cluster_index) {
  MatrixXcd H = H_en(nv, cluster_index);
  if (!nv.no_nn) H += tp(I2,H_nn(nv, cluster_index));
  return H;
}

// nuclear Zeeman Hamiltonian; acts only on cluster
MatrixXcd H_nZ(const nv_system& nv, const uint cluster_index, const Vector3d& gB) {
  const vector<uint> cluster = nv.clusters.at(cluster_index);
  // zero-field splitting and interaction of NV center with magnetic field
  MatrixXcd H = MatrixXcd::Zero(pow(2,cluster.size()),pow(2,cluster.size()));
  for (uint s = 0; s < cluster.size(); s++) {
    // interaction of spin s with magnetic field
    H -= act(dot(gB,I_vec), {s}, cluster.size());
  }
  return H;
}

// perform NV coherence measurement with a static magnetic field
double coherence_measurement(const nv_system& nv, const double w_scan, const double f_DD,
                             const double scan_time) {
  const double w_DD = w_scan/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const vector<double> pulse_times = axy_pulse_times(f_DD,nv.k_DD); // AXY pulse times

  double coherence = 1;
  for (uint cluster = 0; cluster < nv.clusters.size(); cluster++) {
    const uint cluster_size = nv.clusters.at(cluster).size();

    // projections onto |ms> and |0> NV states
    const MatrixXcd proj_m = act(up*up.adjoint(),{0},cluster_size+1); // |ms><ms|
    const MatrixXcd proj_0 = act(dn*dn.adjoint(),{0},cluster_size+1); // |0><0|

    // full system Hamiltonian
    const MatrixXcd H = H_sys(nv,cluster);

    // spin cluster Hamiltonians for single NV states
    const MatrixXcd H_m = ptrace(H*proj_m, {0}); // <ms|H|ms>
    const MatrixXcd H_0 = ptrace(H*proj_0, {0}); // <0|H|0>

    // propagators for sections of the AXY sequence
    const MatrixXcd U1_m = exp(-j*H_m*t_DD*(pulse_times.at(1)-pulse_times.at(0)));
    const MatrixXcd U2_m = exp(-j*H_0*t_DD*(pulse_times.at(2)-pulse_times.at(1)));
    const MatrixXcd U3_m = exp(-j*H_m*t_DD*(pulse_times.at(3)-pulse_times.at(2)));

    const MatrixXcd U1_0 = exp(-j*H_0*t_DD*(pulse_times.at(1)-pulse_times.at(0)));
    const MatrixXcd U2_0 = exp(-j*H_m*t_DD*(pulse_times.at(2)-pulse_times.at(1)));
    const MatrixXcd U3_0 = exp(-j*H_0*t_DD*(pulse_times.at(3)-pulse_times.at(2)));

    // single AXY sequence propagators
    MatrixXcd U_m = U1_m*U2_m*U3_m*U3_0*U2_0*U1_0 * U1_0*U2_0*U3_0*U3_m*U2_m*U1_m;
    MatrixXcd U_0 = U1_0*U2_0*U3_0*U3_m*U2_m*U1_m * U1_m*U2_m*U3_m*U3_0*U2_0*U1_0;

    // propagators for entire scan
    U_m = pow(U_m, int(scan_time/t_DD));
    U_0 = pow(U_0, int(scan_time/t_DD));

    // normalize propagators
    U_m /= sqrt(real(trace(U_m.adjoint()*U_m)) / U_m.rows());
    U_0 /= sqrt(real(trace(U_0.adjoint()*U_0)) / U_0.rows());

    // update coherence
    coherence *= real(trace(U_0.adjoint()*U_m)) / pow(2,cluster_size);
  }
  return coherence;
}

// ---------------------------------------------------------------------------------------
// Control fields and simulation
// ---------------------------------------------------------------------------------------

// return control field for decoupling a single nucleus from other nuclei
control_fields nuclear_decoupling_field(const nv_system& nv, const uint index,
                                        const double phi_rfd, const double theta_rfd) {
  const Vector3d w_j = effective_larmor(nv,index);
  const double w_rfd = w_j.norm()/(1-sin(theta_rfd)/(2*sqrt(2)*nv.scale_factor));
  const double gV_rfd = w_rfd/nv.scale_factor;
  const Vector3d n_rfd = rotate(hat(hyperfine_perp(nv, index)), theta_rfd, hat(w_j));
  return control_fields(gV_rfd*n_rfd, w_rfd, phi_rfd);
}

// perform given rotation on the NV center
protocol rotate_NV(const nv_system& nv, const Vector3d& rotation, const uint spins) {
  if (rotation.squaredNorm() > 0) {
    return protocol(act( exp(-j*dot(rotation,nv.e_S())), {0}, spins), 0, 1);
  } else {
    return protocol::Identity(pow(2,spins));
  }
}

// compute and perform rotation of NV center necessary to generate U_NV
protocol act_NV(const nv_system& nv, const Matrix2cd& U_NV, const uint spins) {
  const Vector4cd H_NV_vec = U_decompose(j*log(U_NV));
  const Vector3d nv_rotation = (xhat*real(H_NV_vec(1))*sqrt(2) +
                                yhat*real(H_NV_vec(2))*nv.ms*sqrt(2) +
                                zhat*real(H_NV_vec(3))*nv.ms*2);
  return rotate_NV(nv,nv_rotation,spins);
}

// simulate propagator with static control fields
protocol simulate_AXY(const nv_system& nv, const uint cluster,
                      const double w_DD, const double f_DD, const axy_harmonic k_DD,
                      const double simulation_time, const double advance_time,
                      const Vector3d gB_ctl) {
  const uint spins = nv.clusters.at(cluster).size()+1;

  // AXY sequence parameters
  const double t_DD = 2*pi/w_DD;
  const double advance = advance_time/t_DD;
  const vector<double> axy_pulses = axy_pulse_times(f_DD,k_DD);
  const vector<double> pulses = advanced_pulse_times(axy_pulses, advance);

  const MatrixXcd H = H_sys(nv, cluster) + H_ctl(nv, cluster, gB_ctl); // full Hamiltonian
  const MatrixXcd X = act_NV(nv, sx, spins).U; // NV center spin flip (pi-)pulse
  MatrixXcd U = MatrixXcd::Identity(H.rows(),H.cols()); // system propagator
  Matrix2cd U_NV = I2; // NV-only propagator
  uint pulse_count = 0;

  // if we need to start with a flipped NV center, flip it
  if (F_AXY(advance, axy_pulses) == -1) {
    U = (X * U).eval();
    U_NV = (sx * U_NV).eval();
  }

  // propagator for whole AXY sequences
  if (simulation_time >= t_DD) {
    MatrixXcd U_AXY = MatrixXcd::Identity(H.rows(),H.cols());
    Matrix2cd U_NV_AXY = I2;
    for (uint i = 1; i < pulses.size(); i++) {
      const double x = pulses.at(i-1);
      const double dx = pulses.at(i) - x;
      U_AXY = (X * exp(-j*dx*t_DD*H) * U_AXY).eval();
      U_NV_AXY = (sx * exp(-j*dx*t_DD*H_NV(nv,gB_ctl)) * U_NV_AXY).eval();
      if(dx > numerical_error) pulse_count++;
    }
    // "undo" the last pulse at x = 1
    U_AXY = (X * U_AXY).eval();
    U_NV_AXY = (sx * U_NV_AXY).eval();
    pulse_count--;

    const long unsigned int sequences = (long unsigned int)(simulation_time/t_DD);
    U = (U_AXY.pow(sequences) * U).eval();
    U_NV = (U_NV_AXY.pow(sequences) * U_NV).eval();
    pulse_count *= sequences;
  }

  // propagator for AXY sequence remainder
  const double remainder = mod(simulation_time/t_DD);
  for (uint i = 1; i < pulses.size(); i++) {
    const double x = pulses.at(i-1);
    const double dx = pulses.at(i) - x;
    if (x + dx < remainder) {
      U = (X * exp(-j*dx*t_DD*H) * U).eval();
      U_NV = (sx * exp(-j*dx*t_DD*H_NV(nv,gB_ctl)) * U_NV).eval();
      if (dx > numerical_error) pulse_count++;
    } else {
      const double dx_f = remainder - x;
      U = (exp(-j*dx_f*t_DD*H) * U).eval();
      U_NV = (exp(-j*dx_f*t_DD*H_NV(nv,gB_ctl)) * U_NV).eval();
      break;
    }
  }
  // rotate into the frame of the NV center and normalize the propagator
  U = (act_NV(nv, U_NV.adjoint(), spins).U * U).eval();
  U /= sqrt(real(trace(U.adjoint()*U)) / U.rows());

  return protocol(U, simulation_time, pulse_count);
}

// simulate propagator with dynamic control fields
protocol simulate_AXY(const nv_system& nv, const uint cluster,
                      const double w_DD, const double f_DD, const axy_harmonic k_DD,
                      const control_fields& controls, const double simulation_time,
                      const double advance_time, const double phi_DD) {
  if (controls.all_fields_static()) {
    return simulate_AXY(nv, cluster, w_DD, f_DD, k_DD, simulation_time,
                        advance_time - phi_DD/w_DD, controls.gB());
  }
  const uint spins = nv.clusters.at(cluster).size()+1;
  const uint D = pow(2,spins);
  if (simulation_time == 0) return protocol::Identity(D);

  // AXY sequence parameters
  const double t_DD = 2*pi/w_DD;
  const double advance = advance_time/t_DD;
  const vector<double> axy_pulses = axy_pulse_times(f_DD,k_DD);
  const vector<double> pulses = advanced_pulse_times(axy_pulses, -phi_DD/(2*pi));

  // largest frequency scale of simulation
  const double frequency_scale = [&]() -> double {
    double largest_control_freq = w_DD;
    Vector3d gB_cap = abs(nv.static_gBz)*zhat;
    for (uint c = 0; c < controls.num(); c++) {
      largest_control_freq = max(largest_control_freq,controls.freqs.at(c));
      gB_cap += (abs(dot(controls.gBs.at(c),xhat))*xhat +
                 abs(dot(controls.gBs.at(c),yhat))*yhat +
                 abs(dot(controls.gBs.at(c),zhat))*zhat);
    }
    return max(largest_control_freq, gB_cap.norm());
  }();


  // integration step number and size
  const uint integration_steps =
    ceil(simulation_time*frequency_scale*nv.integration_factor);
  const double dx = simulation_time/t_DD / integration_steps;

  const MatrixXcd H_0 = H_sys(nv, cluster); // full system Hamiltonian
  const MatrixXcd X = act_NV(nv, sx, spins).U; // NV center spin flip (pi-)pulse
  MatrixXcd U = MatrixXcd::Identity(D,D); // system propagator
  Matrix2cd U_NV = I2; // NV-only propagator
  uint pulse_count = 0;

  // if we need to start with a flipped NV center, flip it
  if (F_AXY(advance - phi_DD/(2*pi), axy_pulses) == -1) {
    U = (X * U).eval();
    U_NV = (sx * U_NV).eval();
  }

  for (uint x_i = 0; x_i < integration_steps; x_i++) {
    const double x = x_i*dx + advance; // time normalized to t_DD

    // determine whether to apply an NV pi-pulse
    const uint pulse = [&]() -> uint {
      const double x_AXY = mod(x);
      for (uint p = 1; p < pulses.size()-1; p++) {
        if (pulses.at(p) < x_AXY + dx) {
          if (pulses.at(p) >= x_AXY) {
            return p;

          }
        } else break;
      }
      return 0;
    }();

    // update propagator
    if (!pulse) {
      const Vector3d gB = controls.gB((x+dx/2)*t_DD);
      const MatrixXcd H = H_0 + H_ctl(nv, cluster, gB);

      U = (exp(-j*dx*t_DD*H) * U).eval();
      U_NV = (exp(-j*dx*t_DD*H_NV(nv,gB)) * U_NV).eval();

    } else { // if (pulse)

      double x_0 = x;
      uint next_pulse = pulse;
      bool even_pulses = true;

      // the following loop handles multiple pulses within a time period of dx
      do {
        double dx_0 = mod(pulses.at(next_pulse) - x_0);

        // fix for clipping bug when dx_0 ~= 1 due to numerical error
        if (dx_0 > 0.5) dx_0 -= 1;

        const Vector3d gB = controls.gB((x_0+dx_0/2)*t_DD);
        const MatrixXcd H = H_0 + H_ctl(nv, cluster, gB);

        U = (X * exp(-j*dx_0*t_DD*H) * U).eval();
        U_NV = (sx * exp(-j*dx_0*t_DD*H_NV(nv,gB)) * U_NV).eval();
        if(abs(dx_0) > numerical_error || even_pulses){
          pulse_count++;
          even_pulses = false;
        } else {
          pulse_count--;
          even_pulses = true;
        }

        x_0 += dx_0;
        next_pulse = next_pulse % (pulses.size()-2) + 1;
      } while (mod(pulses.at(next_pulse) - x) <= dx);

      const double x_f = (x_i+1)*dx + advance;
      const double dx_f = mod(x_f - x_0);
      const Vector3d gB = controls.gB((x_f-dx_f/2)*t_DD);
      const MatrixXcd H = H_0 + H_ctl(nv, cluster, gB);

      U = (exp(-j*dx_f*t_DD*H) * U).eval();
      U_NV = (exp(-j*dx_f*t_DD*H_NV(nv,gB)) * U_NV).eval();
    }
  }
  // rotate into the frame of the NV center and normalize the propagator
  U = (act_NV(nv, U_NV.adjoint(), spins).U * U).eval();
  U /= sqrt(real(trace(U.adjoint()*U)) / U.rows());

  return protocol(U, simulation_time, pulse_count);
}

