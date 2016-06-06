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

  // operations between NV and ST qubits
  static const MatrixXcd SWAP_NVST;

  // spin propagators; Ua corresponds to a Hamiltonian H = h s_a
  static MatrixXcd Ux(const double ht);
  static MatrixXcd Uy(const double ht);
  static MatrixXcd Uz(const double ht);

  // rotation operators
  static MatrixXcd Rx(const double phi);
  static MatrixXcd Ry(const double phi);
  static MatrixXcd Rz(const double phi);

  // spin coupling propagators; Uab corresponds to a Hamiltonian H = h s_a^0 s_b^1
  static MatrixXcd Uxx(const double ht);
  static MatrixXcd Uxy(const double ht);
  static MatrixXcd Uxz(const double ht);
  static MatrixXcd Uyx(const double ht);
  static MatrixXcd Uyy(const double ht);
  static MatrixXcd Uyz(const double ht);
  static MatrixXcd Uzx(const double ht);
  static MatrixXcd Uzy(const double ht);
  static MatrixXcd Uzz(const double ht);

  // controlled phase rotations
  static MatrixXcd cRuu(const double phi);
  static MatrixXcd cRud(const double phi);
  static MatrixXcd cRdu(const double phi);
  static MatrixXcd cRdd(const double phi);
  static MatrixXcd cR(const double phi);

};
