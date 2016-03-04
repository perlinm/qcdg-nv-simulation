#pragma once

using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

//--------------------------------------------------------------------------------------------
// Coordinate systems
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
inline Vector3d natural_axis(const nv_system& nv, const uint index, const double azimuth,
                             const double polar = pi/2){
  return axis(azimuth, polar, natural_basis(nv,index));
}

//--------------------------------------------------------------------------------------------
// General control methods
//--------------------------------------------------------------------------------------------

// polarize an arbitrary state into the pure state psi; warning: not a unitary operation
inline MatrixXcd polarize(const VectorXcd psi){
  return psi*VectorXcd::Ones(psi.size()).adjoint();
}

// propagator U = exp(-i * rotation_angle * sigma_{axis}^{index})
MatrixXcd U_ctl(const nv_system& nv, const uint target, const double target_azimuth,
                const double rotation_angle, const bool exact, const bool adjust_AXY = true,
                const double z_phase = 0);

// compute and perform operationc necessary to act U on target nucleus
MatrixXcd act_target(const nv_system& nv, const uint target, const Matrix2cd& U,
                     const bool exact = false, const bool adjust_AXY = true);

// perform given rotation on a target nucleus
MatrixXcd rotate_target(const nv_system& nv, const uint target, const Vector3d& rotation,
                        const bool exact = false, const bool adjust_AXY = true);

// propagator U = exp(-i * rotation_angle * sigma_{n_1}^{NV}*sigma_{n_2}^{target})
MatrixXcd U_int(const nv_system& nv, const uint target, const Vector3d& nv_axis,
                const double target_azimuth, const double rotation_angle, const bool exact);

//--------------------------------------------------------------------------------------------
// Specific operations
//--------------------------------------------------------------------------------------------

// iSWAP operation
MatrixXcd iSWAP(const nv_system& nv, const uint target, const bool exact);

// SWAP_NVST operation
MatrixXcd SWAP_NVST(const nv_system& nv, const uint idx1, const uint idx2, const bool exact);
