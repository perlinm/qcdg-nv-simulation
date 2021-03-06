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
    const uint subsystem_index = get_index_in_subsystem(nv, index);
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
    rotation = (act(index_rotation, {subsystem_index}, spins) * rotation).eval();
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
  double h_ctl = -gB_ctl/2;
  const double t_rot = abs(2*pi/h_ctl);

  // time for which to apply the control field
  double control_time = mod(angle/h_ctl, t_rot);
  if (control_time > t_rot/2) {
    gB_ctl *= -1;
    h_ctl *= -1;
    control_time = t_rot - control_time;
  }

  const Vector3d axis_ctl = axis(pi/2, target_azimuth, natural_basis(nv,target));
  const control_fields controls(gB_ctl*axis_ctl, w_larmor); // control field object

  protocol U_rotate;
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

    const protocol U_leading = simulate_AXY(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                            controls, leading_time);
    const protocol U_trailing = simulate_AXY(nv, cluster, w_DD_adjusted, f_DD, k_DD,
                                             controls, trailing_time, leading_time);
    U_rotate = U_leading * pow(U_trailing*U_leading, cycles);
  }

  const double flush_time = mod(-control_time - z_angle/w_larmor, t_larmor);
  const protocol U_flush =
    simulate_AXY(nv, cluster, w_DD, f_DD, k_DD, flush_time, control_time);

  return U_flush * U_rotate;
}

// determine and simulate operations necessary to act U on target nucleus
protocol act_target(const nv_system& nv, const uint target, const Matrix2cd& U,
                    const bool exact, const bool adjust_AXY) {

  const uint cluster = get_cluster_containing_target(nv,target);
  const uint spins = nv.clusters.at(cluster).size()+1;
  const uint D = pow(2,spins);

  if (exact) {
    const Matrix2cd R = rotate({xhat,yhat,zhat}, natural_basis(nv,target));
    return protocol(R.adjoint() * U * R);
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

  if (pole_rotation.time < equatorial_rotation.time) return pole_rotation;
  else return equatorial_rotation;
}

// propagator U = exp(-i * angle * sigma_{NV}^{n_1}*I_{target}^{n_2})
protocol U_int(const nv_system& nv, const uint target, const double angle,
               const Vector3d& nv_axis, const double target_azimuth,
               const bool decouple) {
  // identify cluster of target nucleus
  const uint cluster = get_cluster_containing_target(nv,target);
  const uint spins = nv.clusters.at(cluster).size()+1;

  // effective larmor frequency and precession axis,
  //   perpendicular component of hyperfine field,
  //   and maximum possible value of f_DD
  const double w_larmor = effective_larmor(nv,target).norm();
  const Vector3d A_perp = hyperfine_perp(nv,target);
  const double f_DD_max = axy_f_max(nv.k_DD);

  // AXY sequence parameters
  const double w_DD = w_larmor/nv.k_DD; // AXY protocol angular frequency
  const double t_DD = 2*pi/w_DD; // AXY protocol period
  const double dw_min = larmor_resolution(nv,target);
  double f_DD = min(dw_min/(A_perp.norm()*nv.scale_factor), f_DD_max);

  // control field and frequency of precession about its axis
  control_fields controls;
  double w_ctl = 0;

  // default values for AXY sequence phase and effective hyperfine coupling strength
  double phi_DD = target_azimuth;
  double A_int = A_perp.norm();

  // identify index of nucleus to decouple with a control field
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

  // if we have identified a nucleus to decouple
  if (decouple_index >= 0) {
    // precession axis
    const Vector3d w_hat = hat(effective_larmor(nv,target));

    // set control field and precession frequency
    const double gB_ctl =
      min({ exp( (log(2*w_larmor) + log(f_DD_max*A_perp.norm())) / 2),
            2*w_larmor / nv.scale_factor,
            max_gB_ctl});
    const Vector3d axis_ctl = hat(rotate(A_perp, target_azimuth, w_hat));
    controls.add(gB_ctl*axis_ctl, w_larmor);
    w_ctl = gB_ctl/2;

    // set AXY sequence phase
    const Vector3d axis_ctl_rot = rotate(axis_ctl, pi/2, w_hat);
    const Vector3d A_perp_dec = hyperfine_perp(nv, decouple_index);
    const double A_perp_dec_ctl = dot(A_perp_dec, axis_ctl);
    const double A_perp_dec_ctl_rot = dot(A_perp_dec, axis_ctl_rot);
    phi_DD = pi/2 - atan2(A_perp_dec_ctl_rot, A_perp_dec_ctl);

    // set effective hyperfine coupling strength
    const Vector3d A_perp_rot = rotate(A_perp, phi_DD, w_hat);
    A_int = dot(A_perp_rot, axis_ctl);

    // enforce additional limit on allowed values of f_DD
    f_DD = min(f_DD, gB_ctl/(A_perp.norm()*nv.scale_factor));
  }

  // coupling strength and rotation period
  double h_int = f_DD*nv.ms*A_int/4;
  const double t_rot = abs(2*pi/h_int);

  // time for which to interact
  double interaction_time = mod(angle/h_int, t_rot);
  if (interaction_time > t_rot/2) {
    f_DD *= -1;
    h_int *= -1;
    interaction_time = t_rot - interaction_time;
  }

  const unsigned long int cycles = (unsigned long int)(interaction_time/t_DD);
  const double leading_time = interaction_time - cycles*t_DD;

  const protocol U_leading = simulate_AXY(nv, cluster, w_DD, f_DD, nv.k_DD,
                                          controls, leading_time, 0, phi_DD);
  const protocol U_AXY = [&]() -> protocol {
    if (cycles > 0) {
      const double trailing_time = t_DD - leading_time;
      const protocol U_trailing = simulate_AXY(nv, cluster, w_DD, f_DD, nv.k_DD,
                                               controls, trailing_time, leading_time,
                                               phi_DD);
      return U_leading * pow(U_trailing*U_leading, cycles);
    } else return U_leading;
  }();

  // rotate NV coupling axis into its interaction axis (i.e. zhat)
  const protocol nv_axis_rotation = act_NV(nv, rotate(zhat,nv_axis), spins);
  const protocol U_coupling = nv_axis_rotation.adjoint() * U_AXY * nv_axis_rotation;

  // correct for larmor precession of the nucleus
  const double z_angle = mod(interaction_time*w_larmor, 2*pi);
  const double xy_angle = mod(interaction_time*w_ctl, 2*pi);
  const Vector3d xy_axis = rotate(xhat, target_azimuth, zhat);
  const protocol flush_target =
    act_target(nv, target, rotate(xy_angle,xy_axis)*rotate(z_angle,zhat));

  return flush_target * U_coupling;
}

// perform given NV coupling operation on a target nucleus
protocol couple_target(const nv_system& nv, const uint target, const double angle,
                       const Vector3d& nv_axis, const Vector3d& target_axis,
                       const bool exact, const bool decouple, const bool adjust_AXY) {
  if (exact) {
    const Matrix4cd U = exp(-j * angle * tp(dot(s_vec,hat(nv_axis)),
                                            dot(I_vec,hat(target_axis))));
    if (decouple) {
      const Matrix4cd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,target)), {1}, 2);
      return protocol(R.adjoint() * U * R);
    }
    else {
      const uint cluster = get_cluster_containing_target(nv,target);

      vector<uint> nuclei = {};
      for (uint ss = 0; ss < nv.clusters.at(cluster).size(); ss++) {
        const uint nucleus = nv.clusters.at(cluster).at(ss);
        if (is_larmor_pair(nv, target, nucleus)) nuclei.push_back(nucleus);
      }
      const uint spins = nuclei.size()+1;
      assert(spins == 2 || spins == 3);

      MatrixXcd G = MatrixXcd::Identity(pow(2,spins),pow(2,spins));
      for (uint nn = 0; nn < nuclei.size(); nn++) {
        const Matrix4cd R =
          act(rotate({xhat,yhat,zhat}, natural_basis(nv,nuclei.at(nn))), {1}, 2);
        G *= act(R.adjoint() * U * R, {0,nn+1}, spins);
      }
      return protocol(G);
    }
  }

  // compute pitch and azimuth of target axis
  const double target_pitch = real(asin(-j*j*dot(hat(target_axis),zhat)));
  const double target_azimuth = atan2(dot(target_axis,yhat),dot(target_axis,xhat));

  const protocol target_to_equator =
    U_ctl(nv, target, target_pitch, target_azimuth+pi/2, adjust_AXY);
  const protocol coupling_xy =
    U_int(nv, target, angle, nv_axis, target_azimuth, decouple);

  return target_to_equator.adjoint() * coupling_xy * target_to_equator;
}

// ---------------------------------------------------------------------------------------
// Specific operations
// ---------------------------------------------------------------------------------------

// iSWAP operation
protocol iSWAP(const nv_system& nv, const uint target, const bool exact) {
  assert(can_address(nv,target));
  if (exact) {
    const MatrixXcd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,target)), {1}, 2);
    return protocol(R.adjoint() * gates::iSWAP * R);
  } else {
    return (couple_target(nv, target, -pi/2, xhat, xhat) *
            couple_target(nv, target, -pi/2, yhat, yhat));
  }
}

// SWAP operation
protocol SWAP(const nv_system& nv, const uint target, const bool exact) {
  assert(can_address(nv,target));
  if (exact) {
    const MatrixXcd R = act(rotate({xhat,yhat,zhat}, natural_basis(nv,target)), {1}, 2);
    return protocol(R.adjoint() * gates::SWAP * R);
  } else {
    return (couple_target(nv, target, -pi/2, xhat, xhat) *
            couple_target(nv, target, -pi/2, yhat, yhat) *
            couple_target(nv, target, -pi/2, zhat, zhat));
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
    const vector<Vector3d> idx1_basis = natural_basis(nv, idx1);
    const Vector3d x1 = idx1_basis.at(0);
    const Vector3d y1 = idx1_basis.at(1);
    const Vector3d z1 = idx1_basis.at(2);
    const Matrix2cd R1 = rotate({xhat,yhat,zhat}, {x1,y1,z1});
    const Matrix2cd R2 = rotate({xhat,yhat,zhat}, {-y1,x1,z1});
    const MatrixXcd R = act(tp(R1,R2), {1,2}, 3);

    return protocol(R.adjoint() * gates::SWAP_NVST * R);

  } else {
    // compute actual realization of the SWAP_NVST gate
    const protocol hY_NV = rotate_NV(nv, pi/2, yhat, spins);
    const protocol Z_NV = rotate_NV(nv, pi, zhat, spins);
    const protocol XmY_NV = rotate_NV(nv, pi, xhat-yhat, spins);

    const protocol cNOT_AC_NV_adapted =
      hY_NV *
      couple_target(nv, idx1, -pi/2, zhat, xhat) *
      rotate_target(nv, idx1, pi/2, yhat) *
      XmY_NV;

    const Vector3d y1_in_idx2_basis = to_basis(nv, idx2) * from_basis(nv, idx1) * yhat;
    const protocol cNOT_NV_AC_adapted =
      Z_NV *
      rotate_target(nv, idx1, -pi/2, yhat) *
      couple_target(nv, idx1, pi/2, zhat, yhat) *
      couple_target(nv, idx2, pi/2, zhat, y1_in_idx2_basis);

    return (cNOT_AC_NV_adapted.adjoint() *
            cNOT_NV_AC_adapted *
            cNOT_AC_NV_adapted);
  }
}

// identity operation on the cluster containing a target nucleus
//   assumes that the NV center is polarized to the |0> state; acts only on cluster
protocol target_identity(const nv_system& nv, const uint target, const double time,
                         const bool exact, const bool targeting_pair) {
  if (exact) {
    const MatrixXcd H = [&]() -> MatrixXcd {
      const Matrix2cd H_single = -dot(nv.static_gBz*zhat, I_vec);
      if (!targeting_pair) return H_single;
      else return tp(H_single,I2) + tp(I2,H_single);
    }();
    return exp(-j*time*H);
  }

  const uint cluster = get_cluster_containing_target(nv,target);
  MatrixXcd H = H_nZ(nv, cluster, nv.static_gBz*zhat) + H_nn(nv, cluster);

  const MatrixXcd U = exp(-j*time*H);
  return protocol(U, time);
}

// operation to deterministically initialize a thermalized nucleus into |u> or |d>
protocol initialize_spin(const nv_system& nv, const uint target, const bool exact) {
  const uint spins = [&]() -> uint {
    if (exact) return 2;
    else {
      const uint cluster = get_cluster_containing_target(nv,target);
      return nv.clusters.at(cluster).size()+1;
    }
  }();
  return (couple_target(nv, target, pi/2, zhat, yhat, exact) *
          rotate_NV(nv, pi/2, xhat, spins) *
          couple_target(nv, target, pi/2, zhat, xhat, exact) *
          rotate_NV(nv, pi/2, yhat, spins));
}

// operation to probabalistically initialize a thermalized nucleus into |u> +/- |d>
protocol initialize_spin_X(const nv_system& nv, const uint target, const bool exact) {
  const uint spins = [&]() -> uint {
    if (exact) return 2;
    else {
      const uint cluster = get_cluster_containing_target(nv,target);
      return nv.clusters.at(cluster).size()+1;
    }
  }();
  return (rotate_NV(nv, pi/2, xhat, spins) *
          couple_target(nv, target, pi/2, zhat, xhat, exact) *
          rotate_NV(nv, pi/2, yhat, spins));
}

// operation to probabalistically initialize a larmor pair from |dd> into |ud> +/- |du>
protocol initialize_larmor_qubit(const nv_system& nv, const uint idx1, const uint idx2,
                                 const bool exact) {
  assert(is_larmor_pair(nv,idx1,idx2));
  const uint spins = [&]() -> uint {
    if (exact) return 3;
    else {
      const uint cluster = get_cluster_containing_target(nv,idx1);
      assert(in_vector(idx2,nv.clusters.at(cluster)));
      return nv.clusters.at(cluster).size()+1;
    }
  }();
  const bool decouple_spins = false;
  return (rotate_NV(nv, pi/2, yhat, spins) *
          couple_target(nv, idx1, pi/2, zhat, yhat, exact, decouple_spins) *
          rotate_NV(nv, pi/2, yhat, spins));
}
