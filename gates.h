#pragma once
#define EIGEN_USE_MKL_ALL

#include <eigen3/Eigen/Dense> // linear algebra library

using namespace Eigen;

struct gates {

  // single qubit gates
  static const MatrixXcd Z;
  static const MatrixXcd X;
  static const MatrixXcd Y;
  static const MatrixXcd HG;

  // controlled bit-flip and controlled-NOT gates
  static const MatrixXcd cZ;
  static const MatrixXcd cNOT;

  // sqrt(iSWAP), iSWAP, sqrt(SWAP), and SWAP gates
  static const MatrixXcd riSWAP;
  static const MatrixXcd iSWAP;
  static const MatrixXcd rSWAP;
  static const MatrixXcd SWAP;

  // entanglement operator: ud -> S; du -> T
  static const MatrixXcd E;

  // change between coupled and uncoupled basis for two spins
  static const MatrixXcd uncouple;
  static const MatrixXcd couple;

  // singlet/triplet states
  static const VectorXcd S;
  static const VectorXcd T;

  // operations between NV and ST qbits
  static const MatrixXcd SWAP_NVST;

  // rotation operators; Ra corresponds to a Hamiltonian H = h I^a
  static MatrixXcd Rx(const double phi);
  static MatrixXcd Ry(const double phi);
  static MatrixXcd Rz(const double phi);

  // conditional rotation operators; Rab corresponds to a Hamiltonian H = h s_0^a I_1^b
  static MatrixXcd Rxx(const double phi);
  static MatrixXcd Rxy(const double phi);
  static MatrixXcd Rxz(const double phi);
  static MatrixXcd Ryx(const double phi);
  static MatrixXcd Ryy(const double phi);
  static MatrixXcd Ryz(const double phi);
  static MatrixXcd Rzx(const double phi);
  static MatrixXcd Rzy(const double phi);
  static MatrixXcd Rzz(const double phi);

  // controlled phase rotations
  static MatrixXcd cRuu(const double phi);
  static MatrixXcd cRud(const double phi);
  static MatrixXcd cRdu(const double phi);
  static MatrixXcd cRdd(const double phi);
  static MatrixXcd cR(const double phi);

};
