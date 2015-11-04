#pragma once

#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "qp-math.h"
#include "nv-math.h"

// clean up matrix for human readability
inline void clean(MatrixXcd &M){
  remove_artifacts(M);
  remove_phase(M);
  remove_artifacts(M);
}

// returns basis element p for Hamiltonians of a system with N spins
MatrixXcd H_basis_element(int p, int N){
  MatrixXcd spins[4] = {st, sx, sy, sz};
  MatrixXcd b_p = spins[int_bit(p,0)+2*int_bit(p,1)];
  for(int n = 1; n < N; n++){
    b_p = tp(b_p,spins[int_bit(p,2*n)+2*int_bit(p,2*n+1)]);
  }
  return b_p;
}
string H_basis_element_text(int p, int N){
  char spins[4] = {'I','X','Y','Z'};
  stringstream stream;
  for(int n = 0; n < N; n++){
    stream << spins[int_bit(p,2*n)+2*int_bit(p,2*n+1)];
  }
  return stream.str();
}

// flatten matrix into a 1-D vector
inline MatrixXcd flatten(MatrixXcd M){
  M.resize(M.size(),1);
  return M;
}

// returns matrix whose columns are basis Hamiltonians for a system of N spins
MatrixXcd H_basis_matrix(int N){
  MatrixXcd spins[4] = {st, sx, sy, sz};
  MatrixXcd out = MatrixXcd::Zero(pow(4,N),pow(4,N));
  for(int p = 0; p < pow(4,N); p++){
    out.col(p) = flatten(H_basis_element(p,N));
  }
  return out;
}

// decompose Hamiltonian into its basis elements
MatrixXcd H_decompose(const MatrixXcd H, bool fast = true){
  int N = log2(H.rows());
  if(fast) return H_basis_matrix(N).householderQr().solve(flatten(H));
  else return H_basis_matrix(N).fullPivLu().solve(flatten(H));
}

// print Hamiltonian in human-readable form
void H_print(const MatrixXcd H, bool fast = true){
  int N = log2(H.rows());
  MatrixXcd hs = H_decompose(H,fast);
  clean(hs);
  for(int p = 0; p < pow(4,N); p++){
    if(abs(hs(p)) != 0){
      cout << H_basis_element_text(p,N) << ": " << hs(p) << endl;
    }
  }
}

// print state vector in human readable form
void state_print(const MatrixXcd psi){
  int N = psi.size();
  int qbits = log2(N);
  for(int n = 0; n < N; n++){
    if(abs(psi(n)) != 0){
      cout << "|";
      for(int q = 0; q < qbits; q++){
        cout << (qbit_state(q,qbits,n)?"d":"u");
      }
      cout << "> " << psi(n) << endl;
    }
  }
}

// print matrix in human readable form
void matrix_print(const MatrixXcd M){
  int qbits = log2(M.rows());
  for(int m = 0; m < M.rows(); m++){
    for(int n = 0; n < M.cols(); n++){
      if(abs(M(m,n)) != 0){
        cout << "|";
        for(int q = 0; q < qbits; q++){
          cout << (qbit_state(m,qbits,q)?"d":"u");
        }
        cout << "><";
        for(int q = 0; q < qbits; q++){
          cout << (qbit_state(n,qbits,q)?"d":"u");
        }
        cout << "| " << M(m,n) << endl;
      }
    }
  }
}
