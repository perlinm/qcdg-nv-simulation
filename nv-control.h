#pragma once
#define EIGEN_USE_MKL_ALL

#include <eigen3/Eigen/Dense> // linear algebra library

#include "constants.h"
#include "qp-math.h"
#include "gates.h"
#include "nv-math.h"

using namespace std;
using namespace Eigen;

// ---------------------------------------------------------------------------------------
// Coordinate systems
// ---------------------------------------------------------------------------------------

// return "natural" basis of a nucleus
vector<Vector3d> natural_basis(const nv_system& nv, const uint index);

// convert vector from the natural basis of a target into the standard basis
Matrix3d from_basis(const nv_system& nv, const uint target);

// convert vector from the standard basis into the natural basis of a target
inline Matrix3d to_basis(const nv_system& nv, const uint target) {
  return from_basis(nv, target).transpose();
}

// return axis with given azimuth and polar angles in a given basis
inline Vector3d axis(const double polar, const double azimuth,
                     const vector<Vector3d> basis = {xhat, yhat, zhat}) {
  return cos(polar) * basis.at(2) + sin(polar) * ( cos(azimuth) * basis.at(0) +
                                                   sin(azimuth) * basis.at(1) );
}

// rotate into the natural frames of nuclei in a given cluster
// if a target is specified, rotate all nulcei into the frame of that target;
//  otherwise, rotate each nucleus into its own frame
MatrixXcd to_natural_frames(const nv_system& nv, const vector<uint> cluster,
                            const uint target = -1);
inline MatrixXcd to_natural_frames(const nv_system& nv, const uint spin_in_cluster,
                                   const uint target = -1) {
  const uint cluster = get_cluster_containing_target(nv,spin_in_cluster);
  return to_natural_frames(nv, nv.clusters.at(cluster), target);
}

// ---------------------------------------------------------------------------------------
// General control methods
// ---------------------------------------------------------------------------------------

// polarize an arbitrary state into the pure state psi; warning: not a unitary operation
inline MatrixXcd polarize(const VectorXcd psi) {
  return psi*VectorXcd::Ones(psi.size()).adjoint();
}

// propagator U = exp(-i * angle * I_{target}^{axis})
protocol U_ctl(const nv_system& nv, const uint target, const double phase,
               const double target_azimuth, const bool adjust_AXY = true,
               const double z_angle = 0);

// compute and perform operationc necessary to act U on target nucleus
protocol act_target(const nv_system& nv, const uint target, const Matrix2cd& U,
                    const bool exact = false, const bool decouple = false,
                    const bool adjust_AXY = true);

// perform given rotation on a target nucleus
inline protocol rotate_target(const nv_system& nv, const uint target, const double angle,
                              const Vector3d& axis, const bool exact = false,
                              const bool decouple = false, const bool adjust_AXY = true) {
  return act_target(nv, target, rotate(angle, axis), exact, decouple, adjust_AXY);
}

// propagator U = exp(-i * phase * I_{NV}^{n_1}*I_{target}^{n_2})
protocol U_int(const nv_system& nv, const uint target, const double phase,
               const Vector3d& nv_axis, const double target_azimuth,
               const bool decouple = true);

// perform given NV coupling operation on a target nucleus
protocol couple_target(const nv_system& nv, const uint target, const double phase,
                       const Vector3d& nv_axis, const Vector3d& target_axis,
                       const bool exact = false, const bool decouple = true,
                       const bool adjust_AXY = true);

// ---------------------------------------------------------------------------------------
// Specific operations
// ---------------------------------------------------------------------------------------

// iSWAP operation
protocol iSWAP(const nv_system& nv, const uint target, const bool exact = false);

// SWAP operation
protocol SWAP(const nv_system& nv, const uint target, const bool exact = false);

// SWAP operation between NV electron spin and the singlet-triplet (ST) subspace of two
//   nuclear spins; spin bases are: {-z1,y1,x1} for spin 1, and {-y1,x1,z1} for spin 2,
//   where {x1,y1,z1} is the natural basis of spin 1
protocol SWAP_NVST(const nv_system& nv, const uint idx1, const uint idx2,
                   const bool exact = false);
