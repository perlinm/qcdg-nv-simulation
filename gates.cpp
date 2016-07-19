#define EIGEN_USE_MKL_ALL

#include <eigen3/Eigen/Dense> // linear algebra library

#include "constants.h"
#include "qp-math.h"
#include "gates.h"

using namespace Eigen;

// single qubit gates
const MatrixXcd gates::Z = sz;
const MatrixXcd gates::X = sx;
const MatrixXcd gates::Y = sy;
const MatrixXcd gates::HG = (sx+sz)/sqrt(2);

// controlled bit-flip and controlled-NOT gates
const MatrixXcd gates::cZ = (Matrix4cd() <<
                             1, 0, 0, 0,
                             0, 1, 0, 0,
                             0, 0, 1, 0,
                             0, 0, 0,-1).finished();
const MatrixXcd gates::cNOT = (Matrix4cd() <<
                               1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 0, 1,
                               0, 0, 1, 0).finished();

// sqrt(iSWAP), iSWAP, sqrt(SWAP), and SWAP gates
const MatrixXcd gates::riSWAP = (Matrix4cd() <<
                                 1, 0,         0,         0,
                                 0, 1/sqrt(2), j/sqrt(2), 0,
                                 0, j/sqrt(2), 1/sqrt(2), 0,
                                 0, 0,         0,         1).finished();
const MatrixXcd gates::iSWAP = (Matrix4cd() <<
                                1, 0, 0, 0,
                                0, 0, j, 0,
                                0, j, 0, 0,
                                0, 0, 0, 1).finished();
const MatrixXcd gates::rSWAP = (Matrix4cd() <<
                                1, 0,          0,        0,
                                0, (1.+j)/2., (1.-j)/2., 0,
                                0, (1.-j)/2., (1.+j)/2., 0,
                                0, 0,          0,        1).finished();
const MatrixXcd gates::SWAP = (Matrix4cd() <<
                               1, 0, 0, 0,
                               0, 0, 1, 0,
                               0, 1, 0, 0,
                               0, 0, 0, 1).finished();

// entanglement operator: ud -> S; du -> T
const MatrixXcd gates::E = (Matrix4cd() <<
                            1, 0, 0, 1,
                            0, 1, 1, 0,
                            0,-1, 1, 0,
                           -1, 0, 0, 1).finished()/sqrt(2);

// change between coupled and uncoupled basis for two spins
const MatrixXcd gates::uncouple = (Matrix4cd() <<
                                   1, 0,         0,         0,
                                   0, 1/sqrt(2), 1/sqrt(2), 0,
                                   0,-1/sqrt(2), 1/sqrt(2), 0,
                                   0, 0,         0,         1).finished();
const MatrixXcd gates::couple = uncouple.adjoint();

// singlet/triplet states
const VectorXcd gates::S = (ud-du)/sqrt(2);
const VectorXcd gates::T = (ud+du)/sqrt(2);

// operations between NV and ST qbits
const MatrixXcd cNOT_NVST = act(gates::cZ,{0,1},3);
const MatrixXcd cNOT_STNV =
  act(gates::E,{1,2},3) * act(gates::cNOT,{1,0},3) * act(gates::E.adjoint(),{1,2},3);
const MatrixXcd gates::SWAP_NVST = cNOT_NVST * cNOT_STNV * cNOT_NVST;

// rotation operators; Ra corresponds to a Hamiltonian H = h I^a
MatrixXcd gates::Rx(const double phi) { return cos(phi/2)*I2 - j*sin(phi/2)*sx; }
MatrixXcd gates::Ry(const double phi) { return cos(phi/2)*I2 - j*sin(phi/2)*sy; }
MatrixXcd gates::Rz(const double phi) { return cos(phi/2)*I2 - j*sin(phi/2)*sz; }

// conditional rotation operators; Rab corresponds to a Hamiltonian H = h s_0^a I_1^b
MatrixXcd gates::Rxx(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sx,sx); }
MatrixXcd gates::Rxy(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sx,sy); }
MatrixXcd gates::Rxz(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sx,sz); }
MatrixXcd gates::Ryx(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sy,sx); }
MatrixXcd gates::Ryy(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sy,sy); }
MatrixXcd gates::Ryz(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sy,sz); }
MatrixXcd gates::Rzx(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sz,sx); }
MatrixXcd gates::Rzy(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sz,sy); }
MatrixXcd gates::Rzz(const double phi) { return cos(phi/2)*I4 - j*sin(phi/2)*tp(sz,sz); }

// controlled phase rotations
MatrixXcd gates::cRuu(const double phi) {
  return (Matrix4cd() <<
          exp(j*phi), 0, 0, 0,
          0,          1, 0, 0,
          0,          0, 1, 0,
          0,          0, 0, 1).finished();
}
MatrixXcd gates::cRud(const double phi) {
  return (Matrix4cd() <<
          1, 0,          0, 0,
          0, exp(j*phi), 0, 0,
          0, 0,          1, 0,
          0, 0,          0, 1).finished();
}
MatrixXcd gates::cRdu(const double phi) {
  return (Matrix4cd() <<
          1, 0, 0,          0,
          0, 1, 0,          0,
          0, 0, exp(j*phi), 0,
          0, 0, 0,          1).finished();
}
MatrixXcd gates::cRdd(const double phi) {
  return (Matrix4cd() <<
          1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, exp(j*phi)).finished();
}
MatrixXcd gates::cR(const double phi) { return cRdd(phi); }
