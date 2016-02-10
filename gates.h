#pragma once

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "qp-math.h"

struct gates{

  // single qubit gates
  const MatrixXcd Z = sz;
  const MatrixXcd X = sx;
  const MatrixXcd HG = (X+Z)/sqrt(2);

  // controlled bit-flip and controlleg-NOT gates
  const MatrixXcd cZ = (Matrix4cd() << 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1).finished();
  const MatrixXcd cNOT = (Matrix4cd() << 1,0,0,0, 0,1,0,0, 0,0,0,1, 0,0,1,0).finished();

  // sqrt(iSWAP), iSWAP, sqrt(SWAP), and SWAP gates
  const MatrixXcd riSWAP = (Matrix4cd() <<
                            1, 0,         0,         0,
                            0, 1/sqrt(2), j/sqrt(2), 0,
                            0, j/sqrt(2), 1/sqrt(2), 0,
                            0, 0,         0,         1).finished();
  const MatrixXcd iSWAP = (Matrix4cd() << 1,0,0,0, 0,0,j,0, 0,j,0,0, 0,0,0,1).finished();
  const MatrixXcd rSWAP = (Matrix4cd() <<
                           1, 0,         0,         0,
                           0, (1.+j)/2., (1.-j)/2., 0,
                           0, (1.-j)/2., (1.+j)/2., 0,
                           0, 0,         0,         1).finished();
  const MatrixXcd SWAP = (Matrix4cd() << 1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1).finished();

  // entanglement operator: ud -> S; du -> T
  const MatrixXcd E = (Matrix4cd() <<
                       1, 0, 0, 1,
                       0, 1, 1, 0,
                       0,-1, 1, 0,
                       -1, 0, 0, 1).finished()/sqrt(2);
  const MatrixXcd uE = E.inverse();

  // change between coupled and uncoupled basis for two spins
  const MatrixXcd uncouple = (Matrix4cd() <<
                              1, 0,         0,         0,
                              0, 1/sqrt(2), 1/sqrt(2), 0,
                              0,-1/sqrt(2), 1/sqrt(2), 0,
                              0, 0,         0,         1).finished();
  const MatrixXcd couple = uncouple.inverse();

  const VectorXcd S = (ud-du)/sqrt(2);
  const VectorXcd T = (ud+du)/sqrt(2);

  // identity on Hilbert space of two qubits mod the ST subspace
  const MatrixXcd unused_ST_I = uu*uu.transpose() + dd*dd.transpose();

  // phase-flip, bit-flip, and Hadamard gates on an ST qubit
  const MatrixXcd Z_ST = S*S.transpose() - T*T.transpose() + j*unused_ST_I;
  const MatrixXcd X_ST = S*T.transpose() + T*S.transpose() + j*unused_ST_I;
  const MatrixXcd H_ST = ((S+T)*S.transpose()+(S-T)*T.transpose())/sqrt(2.) + j*unused_ST_I;

  // operation to store ST state into NV spin
  const MatrixXcd R_NVST = act(X*HG,{0},3)*act(SWAP,{0,1},3);

  // SWAP between NV and ST qubits
  const MatrixXcd SWAP_NVST = act(E,{1,2},3) * act(cNOT,{1,2},3) * act(X*HG,{0},3)
    * act(cNOT,{0,2},3) * act(SWAP,{0,1},3);

  // spin propagators; Ua corresponds to a Hamiltonian H = h s_a
  MatrixXcd Ux(const double ht){ return cos(ht)*I2 - j*sin(ht)*sx; }
  MatrixXcd Uy(const double ht){ return cos(ht)*I2 - j*sin(ht)*sy; }
  MatrixXcd Uz(const double ht){ return cos(ht)*I2 - j*sin(ht)*sz; }

  // rotation operators
  MatrixXcd Rx(const double phi){ return Ux(phi/2); }
  MatrixXcd Ry(const double phi){ return Uy(phi/2); }
  MatrixXcd Rz(const double phi){ return Uz(phi/2); }

  // spin coupling propagators; Uab corresponds to a Hamiltonian H = h s_a^0 s_b^1
  MatrixXcd Uxx(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sx,sx); }
  MatrixXcd Uxy(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sx,sy); }
  MatrixXcd Uxz(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sx,sz); }
  MatrixXcd Uyx(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sy,sx); }
  MatrixXcd Uyy(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sy,sy); }
  MatrixXcd Uyz(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sy,sz); }
  MatrixXcd Uzx(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sz,sx); }
  MatrixXcd Uzy(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sz,sy); }
  MatrixXcd Uzz(double ht){ return cos(ht)*I4 - j*sin(ht)*tp(sz,sz); }

  // controlled phase rotations
  MatrixXcd cRuu(double phi){
    return (Matrix4cd() << exp(j*phi),0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1).finished();
  }
  MatrixXcd cRud(double phi){
    return (Matrix4cd() << 1,0,0,0, 0,exp(j*phi),0,0, 0,0,1,0, 0,0,0,1).finished();
  }
  MatrixXcd cRdu(double phi){
    return (Matrix4cd() << 1,0,0,0, 0,1,0,0, 0,0,exp(j*phi),0, 0,0,0,1).finished();
  }
  MatrixXcd cRdd(double phi){
    return (Matrix4cd() << 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,exp(j*phi)).finished();
  }
  MatrixXcd cR(double phi){ return cRdd(phi); }

  // rotation operators
  MatrixXcd Rx_ST(double phi){
    return (cos(phi/2.)*(S*S.transpose()+T*T.transpose())
            - j*sin(phi/2.)*(S*T.transpose()+T*S.transpose()) + unused_ST_I);
  }
  MatrixXcd Ry_ST(double phi){
    return (cos(phi/2)*(S*S.transpose()+T*T.transpose())
            + sin(phi/2)*(-S*T.transpose()+T*S.transpose()) + unused_ST_I);
  }
  MatrixXcd Rz_ST(double phi){
    return (exp(-j*phi/2.)*S*S.transpose() +
            exp( j*phi/2.)*T*T.transpose() + unused_ST_I);
  }

};
