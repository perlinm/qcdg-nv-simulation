#pragma once

#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
using namespace Eigen;

#include "qp-math.h"
#include "nv-math.h"

// clean up matrix for human readability
inline MatrixXcd clean(const MatrixXcd& M){
  return remove_artifacts(remove_phase(remove_artifacts(M)));
}

// returns element p of a basis for operators acting on a system with N qubits
string U_basis_element_text(const uint p, const uint N){
  const char spins[4] = {'I','X','Y','Z'};
  stringstream stream;
  for(uint n = 0; n < N; n++){
    stream << spins[int_bit(p,2*n)+2*int_bit(p,2*n+1)];
  }
  return stream.str();
}

// print operator in human-readable form
void U_print(const MatrixXcd& U, const bool fast = true){
  const int N = log2(U.rows());
  const MatrixXcd hs = clean(U_decompose(U,fast));
  for(int p = 0; p < pow(4,N); p++){
    if(abs(hs(p)) != 0){
      cout << U_basis_element_text(p,N) << ": " << hs(p) << endl;
    }
  }
}

// print state vector in human readable form
void state_print(const MatrixXcd& psi){
  const uint N = psi.size();
  const uint qbits = log2(N);
  for(uint n = 0; n < N; n++){
    if(abs(psi(n)) != 0){
      cout << "|";
      for(uint q = 0; q < qbits; q++){
        cout << (qbit_state(q,qbits,n)?"d":"u");
      }
      cout << "> " << psi(n) << endl;
    }
  }
}

// print matrix in human readable form
void matrix_print(const MatrixXcd& M){
  const int qbits = log2(M.rows());
  for(uint m = 0; m < M.rows(); m++){
    for(uint n = 0; n < M.cols(); n++){
      if(abs(M(m,n)) != 0){
        cout << "|";
        for(uint q = 0; q < qbits; q++){
          cout << (qbit_state(m,qbits,q)?"d":"u");
        }
        cout << "><";
        for(uint q = 0; q < qbits; q++){
          cout << (qbit_state(n,qbits,q)?"d":"u");
        }
        cout << "| " << M(m,n) << endl;
      }
    }
  }
}
