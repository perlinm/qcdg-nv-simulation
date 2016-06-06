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

// operations between NV and ST qubits
const MatrixXcd cNOT_NVST = act(gates::cZ,{0,1},3);
const MatrixXcd cNOT_STNV =
  act(gates::E,{1,2},3) * act(gates::cNOT,{1,0},3) * act(gates::E.adjoint(),{1,2},3);
const MatrixXcd gates::SWAP_NVST = cNOT_NVST * cNOT_STNV * cNOT_NVST;

// spin propagators; Ua corresponds to a Hamiltonian H = h s_a
MatrixXcd gates::Ux(const double ht) { return cos(ht)*I2 - j*sin(ht)*sx; }
MatrixXcd gates::Uy(const double ht) { return cos(ht)*I2 - j*sin(ht)*sy; }
MatrixXcd gates::Uz(const double ht) { return cos(ht)*I2 - j*sin(ht)*sz; }

// rotation operators
MatrixXcd gates::Rx(const double phi) { return Ux(phi/2); }
MatrixXcd gates::Ry(const double phi) { return Uy(phi/2); }
MatrixXcd gates::Rz(const double phi) { return Uz(phi/2); }

// spin coupling propagators; Uab corresponds to a Hamiltonian H = h s_a^0 s_b^1
MatrixXcd gates::Uxx(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sx,sx); }
MatrixXcd gates::Uxy(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sx,sy); }
MatrixXcd gates::Uxz(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sx,sz); }
MatrixXcd gates::Uyx(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sy,sx); }
MatrixXcd gates::Uyy(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sy,sy); }
MatrixXcd gates::Uyz(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sy,sz); }
MatrixXcd gates::Uzx(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sz,sx); }
MatrixXcd gates::Uzy(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sz,sy); }
MatrixXcd gates::Uzz(const double ht) { return cos(ht)*I4 - j*sin(ht)*tp(sz,sz); }

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
