#pragma once

using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

// ---------------------------------------------------------------------------------------
// Coordinate systems
// ---------------------------------------------------------------------------------------

// return "natural" basis of a nucleus
vector<Vector3d> natural_basis(const nv_system& nv, const uint index);

// return axis with given azimuth and polar angles in a given basis
inline Vector3d axis(const double polar, const double azimuth,
                     const vector<Vector3d> basis = {xhat, yhat, zhat}){
  return cos(polar) * basis.at(2) + sin(polar) * ( cos(azimuth) * basis.at(0) +
                                                   sin(azimuth) * basis.at(1) );
}

// rotate into the natural frames of all nuclei in the cluster
MatrixXcd to_natural_frames(const nv_system& nv, const vector<uint> cluster);
inline MatrixXcd to_natural_frames(const nv_system& nv, const uint cluster){
  return to_natural_frames(nv, nv.clusters.at(cluster));
}

// ---------------------------------------------------------------------------------------
// General control methods
// ---------------------------------------------------------------------------------------

// polarize an arbitrary state into the pure state psi; warning: not a unitary operation
inline MatrixXcd polarize(const VectorXcd psi){
  return psi*VectorXcd::Ones(psi.size()).adjoint();
}

// propagator U = exp(-i * phase * sigma_{axis}^{index})
protocol U_ctl(const nv_system& nv, const uint target, const double phase,
               const double target_azimuth, const bool adjust_AXY = true,
               const double z_phase = 0);

// compute and perform operationc necessary to act U on target nucleus
protocol act_target(const nv_system& nv, const uint target, const Matrix2cd& U,
                    const bool exact = false, const bool adjust_AXY = true);

// perform given rotation on a target nucleus
inline protocol rotate_target(const nv_system& nv, const uint target,
                              const Vector3d& rotation, const bool exact = false,
                              const bool adjust_AXY = true){
  return act_target(nv, target, rotate(rotation), exact, adjust_AXY);
}

// propagator U = exp(-i * phase * sigma_{n_1}^{NV}*sigma_{n_2}^{target})
protocol U_int(const nv_system& nv, const uint target, const double phase,
               const Vector3d& nv_axis, const double target_azimuth);

// perform given NV coupling operation on a target nucleus
protocol couple_target(const nv_system& nv, const uint target, const double phase,
                       const Vector3d& nv_axis, const Vector3d& target_axis,
                       const bool exact = false, const bool adjust_AXY = true);

// ---------------------------------------------------------------------------------------
// Specific operations
// ---------------------------------------------------------------------------------------

// iSWAP operation
inline protocol iSWAP(const nv_system& nv, const uint index, const bool exact = false){
  return (couple_target(nv, index, -pi/4, xhat, xhat, exact) *
          couple_target(nv, index, -pi/4, yhat, yhat, exact));
}

// SWAP operation
inline protocol SWAP(const nv_system& nv, const uint index, const bool exact = false){
  return (couple_target(nv, index, -pi/4, xhat, xhat, exact) *
          couple_target(nv, index, -pi/4, yhat, yhat, exact) *
          couple_target(nv, index, -pi/4, zhat, zhat, exact));
}

// SWAP_NVST operation
protocol SWAP_NVST(const nv_system& nv, const uint idx1, const uint idx2,
                   const bool exact = false);
