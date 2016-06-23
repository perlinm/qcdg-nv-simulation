#define EIGEN_USE_MKL_ALL

#include <iostream> // for standard output
#include <eigen3/Eigen/Dense> // linear algebra library

#include "constants.h"
#include "qp-math.h"
#include "gates.h"
#include "nv-math.h"
#include "nv-control.h"

using namespace std;
using namespace Eigen;

// ---------------------------------------------------------------------------------------
// Coordinate systems
// ---------------------------------------------------------------------------------------

// return "natural" basis of a nucleus
vector<Vector3d> natural_basis(const nv_system& nv, const uint index) {
  const Vector3d target_zhat = hat(effective_larmor(nv,index));
  const Vector3d target_xhat = hat(hyperfine_perp(nv,index));
  const Vector3d target_yhat = target_zhat.cross(target_xhat);
  return {target_xhat, target_yhat, target_zhat};
}

// convert vector from the standard basis into the natural basis of a target
Matrix3d from_basis(const nv_system& nv, const uint target){
  const vector<Vector3d> target_basis = natural_basis(nv, target);
  return (target_basis.at(0) * xhat.transpose() +
          target_basis.at(1) * yhat.transpose() +
          target_basis.at(2) * zhat.transpose());
}

// rotate into the natural frames of nuclei in a given cluster
// if a target is specified, rotate all nulcei into the frame of that target;
//  otherwise, rotate each nucleus into its own frame
MatrixXcd to_natural_frames(const nv_system& nv, const vector<uint> cluster,
                            const uint target) {
  const uint spins = cluster.size()+1;
  MatrixXcd rotation = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
  for (uint index: cluster) {
    const uint index_in_cluster = get_index_in_cluster(index, cluster);
    const Matrix2cd index_rotation = [&]() -> Matrix2cd {
      if (can_address(nv,index)) {
        if (target == (uint)(-1)) {
          return rotate(natural_basis(nv,index), {xhat,yhat,zhat});
        } else {
          return rotate(natural_basis(nv,target), {xhat,yhat,zhat});
        }
      } else {
        return I2;
      }
    }();
    rotation = (act(index_rotation, {index_in_cluster+1}, spins) * rotation).eval();
  }
  return rotation;
}

// ---------------------------------------------------------------------------------------
// General control methods
// ---------------------------------------------------------------------------------------

// propagator U = exp(-i * angle * I_{target}^{axis})
protocol U_ctl(const nv_system& nv, const uint target, const double angle,
               const double target_azimuth, const bool adjust_AXY, const double z_angle) {
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_target(nv,target);

  // larmor frequency of target nucleus
  const double w_larmor = effective_larmor(nv,target).norm();
  const double t_larmor = 2*pi/w_larmor;

  // AXY protocol parameters
  const double sA = nv.scale_factor * hyperfine(nv,target).norm();
  const double w_DD = [&]() -> double {
    const double w_DD_large = (w_larmor+sA)/3.;
    if (w_larmor < sA) {
      return w_DD_large;
    } else {
      const uint k_m = 2*int(0.5 * (w_larmor/sA-1) );
      const double w_DD_small = (w_larmor-sA)/k_m;
      if (w_DD_small > sA && isfinite(w_DD_small)) {
        return w_DD_small;
      } else {
        return w_DD_large;
      }
    }
  }();

  const axy_harmonic k_DD = abs(w_DD - w_larmor) < abs(3*w_DD - w_larmor) ? first : third;
  const double f_DD = 0;

  const double dw_min = larmor_resolution(nv,target);
  double gB_ctl = min(dw_min/nv.scale_factor, max_gB_ctl);

  // coupling strength and rotation period
  const double h_ctl = -gB_ctl/2;
  const double t_rot = abs(2*pi/h_ctl);

  // time for which to apply the control field
  double control_time = mod(angle/h_ctl, t_rot); // control operation time
  if (control_time > t_rot/2) {
    gB_ctl *= -1;
    control_time = t_rot - control_time;
  }

  const Vector3d axis_ctl = axis(pi/2, target_azimuth, natural_basis(nv,target));
  const control_fields controls(gB_ctl*axis_ctl, w_larmor); // control field object

  MatrixXcd U_rotate;
  if (!adjust_AXY) {
    U_rotate = simulate_AXY(nv, cluster, w_DD, f_DD, k_DD, controls, control_time);
  } else {

    const double freq_ratio = [&]() -> double {
      if (w_DD < w_larmor) {
        return 2*round(0.5*w_larmor/w_DD);
      } else {
        return 1/round(w_DD/w_larmor);
      }
    }();
    const double w_DD_adjusted = w_larmor/freq_ratio;
    const double t_DD_adjusted = 2*pi/w_DD_adjusted;

    const double cycle_time = max(t_DD_adjusted, t_larmor);
    const unsigned long int cycles = (unsigned long int)(control_time/cycle_time);

    const double leading_time = control_time - cycles*cycle_time;
    const double trailing_time = [&]() -> double {
      if (cycles == 0) return 0;
      else return cycle_time - leading_time;
    }();

    const MatrixXcd U_leading = simulate_AXY(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                             controls, leading_time);
    const MatrixXcd U_trailing = simulate_AXY(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                              controls, trailing_time, leading_time);
    U_rotate = U_leading * pow(U_trailing*U_leading, cycles);
  }

  const double flush_time = mod(-control_time - z_angle/w_larmor, t_larmor);
  const MatrixXcd U_flush =
    simulate_AXY(nv, cluster, w_DD, f_DD, k_DD, flush_time, control_time);
  return protocol(U_flush * U_rotate, control_time + flush_time);
}

// determine and simulate operations necessary to act U on target nucleus
protocol act_target(const nv_system& nv, const uint target, const Matrix2cd& U,
                    const bool exact, const bool decouple, const bool adjust_AXY) {
  // when performing this operation on larmor sets, we cannot "decouple",
  //   i.e. act on only one spin
  assert(exact || !decouple);

  const uint cluster = get_cluster_containing_target(nv,target);
  const uint target_in_cluster = get_index_in_cluster(nv,target);
  const uint spins = nv.clusters.at(cluster).size()+1;
  const uint D = pow(2,spins);

  if (exact) {
    // compute rotation to act on target nucleus
    const Matrix2cd to_target_basis = rotate(natural_basis(nv,target), {xhat,yhat,zhat});
    const Matrix2cd rotation = to_target_basis * U * to_target_basis.adjoint();

    if (decouple) {
      // rotate only target nucleus
      return protocol(act(rotation, {target_in_cluster+1}, spins), 0);
    } else {
      // rotate all nuclei with the same effective larmor frequency as the target
      MatrixXcd G = MatrixXcd::Identity(D,D);
      for (uint s = 0; s < nv.clusters.at(cluster).size(); s++) {
        if (is_larmor_pair(nv,target,nv.clusters.at(cluster).at(s))) {
          G *= act(rotation, {s+1}, spins);
        }
      }
      return protocol(G,0);
    }
  }

  const Vector4cd H_vec = U_decompose(j*log(U)*2);
  const double rx = real(H_vec(1));
  const double ry = real(H_vec(2));
  const double rz = real(H_vec(3));

  const double angle = sqrt(rx*rx + ry*ry + rz*rz);
  if (angle < numerical_error || 2*pi - angle < numerical_error) {
    return protocol::Identity(D);
  }

  const double azimuth = atan2(ry,rx);
  const double pitch = asin(rz/angle);

  const protocol pole_rotation = [&]() -> protocol {
    const int pole = pitch > 0 ? 1 : -1; // "north" vs "south" pole
    const double angle_to_pole = pi/2 - abs(pitch);

    const protocol to_pole = [&]() -> protocol {
      if (angle_to_pole < numerical_error) return protocol::Identity(D);
      else return U_ctl(nv, target, pole*angle_to_pole, azimuth-pi/2, adjust_AXY);
    }();
    const protocol rotate_z = U_ctl(nv, target, 0, 0, adjust_AXY, pole*angle);

    return to_pole.adjoint() * rotate_z * to_pole;
  }();

  const protocol equatorial_rotation = [&]() -> protocol {
    const protocol to_equator = [&]() -> protocol {
      if (abs(pitch) < numerical_error) return protocol::Identity(D);
      else return U_ctl(nv, target, pitch, azimuth+pi/2, adjust_AXY);
    }();
    const protocol rotate_xy = U_ctl(nv, target, angle, azimuth, adjust_AXY);

    return to_equator.adjoint() * rotate_xy * to_equator;
  }();

  if (pole_rotation.t < equatorial_rotation.t) return pole_rotation;
  else return equatorial_rotation;
}

// propagator U = exp(-i * phase * I_{NV}^{n_1}*I_{target}^{n_2})
protocol U_int(const nv_system& nv, const uint target, const double phase,
               const Vector3d& nv_axis, const double target_azimuth,
               const bool decouple) {
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_target(nv,target);
  const uint spins = nv.clusters.at(cluster).size()+1;

  // effective larmor frequency and rotation axis,
  //   perpendicular component of hyperfine field,
  //   and maximum possible value of f_DD
  const double w_larmor = effective_larmor(nv,target).norm();
  const Vector3d w_hat = hat(effective_larmor(nv,target));
  const Vector3d A_perp = hyperfine_perp(nv,target);
  const double f_DD_max = axy_f_max(nv.k_DD);

  // control fields, AXY sequence phase, and effective hyperfine coupling strength
  control_fields controls;
  double phi_DD;
  double A_int;

  const int decouple_index = [&]() -> int {
    if (!decouple) return -1;
    else {
      for (uint index: nv.clusters.at(cluster)) {
        if (index == target) continue;
        if (is_larmor_pair(nv,index,target)) {
          return index;
        }
      }
      return -1;
    }
  }();

  if (decouple_index >= 0) {
    // set control field
    const double gB_ctl =
      min({ exp( (log(2*w_larmor) + log(f_DD_max*A_perp.norm())) / 2),
            2*w_larmor / nv.scale_factor,
            max_gB_ctl});
    const Vector3d axis_ctl = hat(rotate(A_perp, target_azimuth, w_hat));
    controls.add(gB_ctl*axis_ctl, w_larmor);

    // set phi_DD
    const Vector3d axis_ctl_rot = rotate(axis_ctl, pi/2, w_hat);
    const Vector3d A_perp_dec = hyperfine_perp(nv, decouple_index);
    const double A_perp_dec_ctl = dot(A_perp_dec, axis_ctl);
    const double A_perp_dec_ctl_rot = dot(A_perp_dec, axis_ctl_rot);
    phi_DD = pi/2 - atan2(A_perp_dec_ctl_rot, A_perp_dec_ctl);

    // set A_int
    const Vector3d A_perp_rot = rotate(A_perp, phi_DD, w_hat);
    A_int = dot(A_perp_rot, axis_ctl);

  } else { // there is nothing to decouple
    phi_DD = target_azimuth;
    A_int = A_perp.norm();
  }

  // AXY sequence parameters
  const double w_DD = w_larmor/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  double f_DD = min(larmor_resolution(nv,target)/(A_int*nv.scale_factor), f_DD_max);
  if (decouple_index >= 0 &&
      controls.gB().norm() / (f_DD*A_int) < nv.scale_factor) {
    f_DD = controls.gB().norm() / (nv.scale_factor*A_int);
  }

  // coupling strength and rotation period
  const double h_int = f_DD*nv.ms*A_int/2;
  const double t_rot = abs(4*pi/h_int);

  // time for which to interact
  double interaction_time = mod(phase/h_int, t_rot);
  if (interaction_time > t_rot/2) {
    f_DD *= -1;
    interaction_time = t_rot - interaction_time;
  }

  const unsigned long int cycles = (unsigned long int)(interaction_time/t_DD);
  const double leading_time = interaction_time - cycles*t_DD;

  const MatrixXcd U_leading = simulate_AXY(nv, cluster, w_DD, f_DD, nv.k_DD,
                                           controls, leading_time, 0, phi_DD);
  const MatrixXcd U_AXY = [&]() -> MatrixXcd {
    if (cycles > 0) {
      const double trailing_time = t_DD - leading_time;
      const MatrixXcd U_trailing = simulate_AXY(nv, cluster, w_DD, f_DD, nv.k_DD,
                                                controls, trailing_time, leading_time,
                                                phi_DD);
      return U_leading * pow(U_trailing*U_leading, cycles);
    } else return U_leading;
  }();

  // rotate NV coupling axis into its interaction axis (i.e. zhat)
  const MatrixXcd nv_axis_rotation = act_NV(nv, rotate(zhat,nv_axis), spins);
  const MatrixXcd U_coupling = nv_axis_rotation.adjoint() * U_AXY * nv_axis_rotation;

  // correct for larmor precession of the nucleus
  const double z_angle = mod(interaction_time*w_larmor, 2*pi);
  const double w_ctl = controls.gB().norm()/2;
  const double xy_angle = mod(interaction_time*w_ctl, 2*pi);
  const Vector3d xy_axis = rotate(xhat, target_azimuth, zhat);
  const protocol flush_target =
    act_target(nv, target, rotate(xy_angle,xy_axis)*rotate(z_angle,zhat));

  return flush_target * protocol(U_coupling, interaction_time);
}

// perform given NV coupling operation on a target nucleus
protocol couple_target(const nv_system& nv, const uint target, const double phase,
                       const Vector3d& nv_axis, const Vector3d& target_axis,
                       const bool exact, const bool decouple, const bool adjust_AXY) {
  if (exact) {
    // identify cluster of target nucleus
    const uint cluster = get_cluster_containing_target(nv,target);
    const uint target_in_cluster = get_index_in_cluster(nv,target);
    const uint spins = nv.clusters.at(cluster).size()+1;

    const Matrix4cd U = exp(-j * phase * tp(dot(s_vec/2,hat(nv_axis)),
                                            dot(s_vec/2,hat(target_axis))));
    if (decouple) {
      const Matrix4cd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,target)), {1}, 2);
      return protocol(act(R.adjoint() * U * R, {0,target_in_cluster+1}, spins), 0);
    } else {
      MatrixXcd G = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
      for (uint s = 0; s < nv.clusters.at(cluster).size(); s++) {
        const uint n = nv.clusters.at(cluster).at(s);
        if (is_larmor_pair(nv, target, n)) {
          const Matrix4cd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,n)), {1}, 2);
          G *= act(R.adjoint() * U * R, {0,s+1}, spins);
        }
      }
      return protocol(G,0);
    }
  }

  // compute pitch and azimuth of target axis
  const double target_pitch = real(asin(-j*j*dot(hat(target_axis),zhat)));
  const double target_azimuth = atan2(dot(target_axis,yhat),dot(target_axis,xhat));

  const protocol target_to_equator =
    U_ctl(nv, target, target_pitch, target_azimuth+pi/2, adjust_AXY);
  const protocol coupling_xy =
    U_int(nv, target, phase, nv_axis, target_azimuth, decouple);

  return target_to_equator.adjoint() * coupling_xy * target_to_equator;
}

// ---------------------------------------------------------------------------------------
// Specific operations
// ---------------------------------------------------------------------------------------

// iSWAP operation
protocol iSWAP(const nv_system& nv, const uint target, const bool exact) {
  assert(can_address(nv,target));
  if (exact) {
    const uint cluster = get_cluster_containing_target(nv,target);
    const uint target_in_cluster = get_index_in_cluster(nv,target);
    const uint spins = nv.clusters.at(cluster).size()+1;

    const MatrixXcd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,target)), {1}, 2);
    return protocol(act(R.adjoint() * gates::iSWAP * R,
                        {0,target_in_cluster+1}, spins));
  } else {
    return (couple_target(nv, target, -pi, xhat, xhat) *
            couple_target(nv, target, -pi, yhat, yhat));
  }
}

// SWAP operation
protocol SWAP(const nv_system& nv, const uint target, const bool exact) {
  assert(can_address(nv,target));
  if (exact) {
    const uint cluster = get_cluster_containing_target(nv,target);
    const uint target_in_cluster = get_index_in_cluster(nv,target);
    const uint spins = nv.clusters.at(cluster).size()+1;

    const MatrixXcd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,target)), {1}, 2);
    return protocol(act(R.adjoint() * gates::SWAP * R,
                        {0,target_in_cluster+1}, spins));
  } else {
    return (couple_target(nv, target, -pi, xhat, xhat) *
            couple_target(nv, target, -pi, yhat, yhat) *
            couple_target(nv, target, -pi, zhat, zhat));
  }
}

// SWAP operation between NV electron spin and the singlet-triplet (ST) subspace of two
//   nuclear spins; spin 1 is in its natural basis {x1,y1,z1}, whereas spin 2 is in the
//   basis {-y1,x1,z1}
protocol SWAP_NVST(const nv_system& nv, const uint idx1, const uint idx2,
                   const bool exact) {
  // assert that both target nuclei are larmor pairs in the same cluster
  assert(is_larmor_pair(nv,idx1,idx2));
  const uint cluster = get_cluster_containing_target(nv,idx1);
  const uint spins = nv.clusters.at(cluster).size()+1;
  assert(in_vector(idx2,nv.clusters.at(cluster)));

  if (exact) {
    // return exact SWAP_NVST gate in the appropriate bases
    const uint cidx1 = get_index_in_cluster(nv, idx1)+1;
    const uint cidx2 = get_index_in_cluster(nv, idx2)+1;

    const vector<Vector3d> idx1_basis = natural_basis(nv, idx1);
    const Vector3d x1 = idx1_basis.at(0);
    const Vector3d y1 = idx1_basis.at(1);
    const Vector3d z1 = idx1_basis.at(2);
    const Matrix2cd R1 = rotate({xhat,yhat,zhat}, {x1,y1,z1});
    const Matrix2cd R2 = rotate({xhat,yhat,zhat}, {-y1,x1,z1});
    const MatrixXcd R = act(tp(R1,R2), {1,2}, 3);

    return protocol(act(R.adjoint() * gates::SWAP_NVST * R, {0,cidx1,cidx2}, spins));

  } else {
    // compute actual realization of the SWAP_NVST gate
    const protocol qY_NV = protocol(act_NV(nv, rotate(pi/2,yhat), spins));
    const protocol hZ_NV = protocol(act_NV(nv, rotate(pi,zhat), spins));
    const protocol hXmY_NV = protocol(act_NV(nv, rotate(pi,xhat-yhat), spins));

    const protocol cNOT_AC_NV_adapted =
      qY_NV *
      couple_target(nv, idx1, -pi, zhat, xhat) *
      rotate_target(nv, idx1, pi/2, yhat) *
      hXmY_NV;

    const Vector3d y1_in_idx2_basis = to_basis(nv, idx2) * from_basis(nv, idx1) * yhat;
    const protocol cNOT_NV_AC_adapted =
      hZ_NV *
      rotate_target(nv, idx1, -pi/2, yhat) *
      couple_target(nv, idx1, pi, zhat, yhat) *
      couple_target(nv, idx2, pi, zhat, y1_in_idx2_basis);

    return (cNOT_AC_NV_adapted.adjoint() *
            cNOT_NV_AC_adapted *
            cNOT_AC_NV_adapted);
  }
}
