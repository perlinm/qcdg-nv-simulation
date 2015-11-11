using namespace std;

#include <eigen3/Eigen/Dense> // linear algebra library
#include <eigen3/unsupported/Eigen/KroneckerProduct> // provides tensor product
#include <eigen3/unsupported/Eigen/MatrixFunctions> // provides matrix functions
using namespace Eigen;

#include "constants.h"
#include "qp-math.h"

// tensor product of many matrices
MatrixXcd tp(const initializer_list<MatrixXcd> list){
  MatrixXcd out = I1;
  for(MatrixXcd elem: list){
    out = tp(out,elem);
  }
  return out;
}

// remove numerical artifacts from a matrix
void remove_artifacts(MatrixXcd& A, double threshold){
  for(uint m = 0; m < A.rows(); m++){
    for(uint n = 0; n < A.cols(); n++){
      if(abs(A(m,n).real()) < threshold) A(m,n) -= A(m,n).real();
      if(abs(A(m,n).imag()) < threshold) A(m,n) -= A(m,n).imag()*j;
    }
  }
}

// get global phase of matrix
complex<double> get_phase(const MatrixXcd& A){
  for(uint m = 0; m < A.rows(); m++){
    for(uint n = 0; n < A.cols(); n++){
      if(abs(A(m,n)) != 0){
        complex<double> phase = A(m,n)/abs(A(m,n));
        if(phase.real() < 0) phase = -phase;
        return phase;
      }
    }
  }
  return 1;
}

//--------------------------------------------------------------------------------------------
// Operator rearrangement
//--------------------------------------------------------------------------------------------

// generate matrix B to act A on qbits qs_act out of qbits_new
MatrixXcd act(const MatrixXcd& A, const vector<uint> qs_act, uint qbits_new){
  assert(A.rows() == A.cols()); // A should be square

  // number of qbits A acted on
  const uint qbits_old = qs_act.size();
  assert(qbits_old == log2(A.rows()));

  // vector of qubits we are ignoring
  vector<uint> qs_ignore;
  for(uint i = 0; i < qbits_new; i++){
    if(!in_vector(i,qs_act)){
      qs_ignore.push_back(i);
    }
  }

  // initialize B (output) to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries A(m,n)
  for(uint m = 0; m < A.rows(); m++){
    for(uint n = 0; n < A.cols(); n++){

      // get contribution of substates |m-><n-| to indices of B
      uint b_m = 0, b_n = 0;
      for(uint q = 0; q < qbits_old; q++){
        if(qbit_state(q,qbits_old,m)) b_m += bit_int(qs_act.at(q),qbits_new);
        if(qbit_state(q,qbits_old,n)) b_n += bit_int(qs_act.at(q),qbits_new);
      }

      // loop over all elements of the form |ms><ns| in B
      for(uint s = 0; s < pow(2,qs_ignore.size()); s++){
        uint b_out = b_m, b_in = b_n;
        for(uint q = 0; q < qs_ignore.size(); q++){
          if(qbit_state(q,qs_ignore.size(),s)){
            b_out += bit_int(qs_ignore.at(q),qbits_new);
            b_in += bit_int(qs_ignore.at(q),qbits_new);
          }
        }
        B(b_out,b_in) = A(m,n);
      }
    }
  }
  return B;
}

// perform a partial trace over qbits qs_trace
MatrixXcd ptrace(const MatrixXcd& A, const vector<uint> qs_trace){
  assert(A.rows() == A.cols()); // A should be square

  // number of qbits A acted on
  const uint qbits_old = log2(A.rows());
  const uint qbits_new = qbits_old - qs_trace.size();

  // vector of qubits we are keeping
  vector<uint> qs_keep;
  for(uint i = 0; i < qbits_old; i++){
    if(!in_vector(i,qs_trace)) qs_keep.push_back(i);
  }
  assert(qbits_new == qs_keep.size());

  // initialize B (output) to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries B(m,n)
  for(uint m = 0; m < B.rows(); m++){
    for(uint n = 0; n < B.cols(); n++){

      // get contribution of substates |m-><n-| to indices of A
      uint a_m = 0, a_n = 0;
      for(uint q = 0; q < qbits_new; q++){
        if(qbit_state(q,qbits_new,m)) a_m += bit_int(qs_keep.at(q),qbits_old);
        if(qbit_state(q,qbits_new,n)) a_n += bit_int(qs_keep.at(q),qbits_old);
      }

      // loop over all elements of the form |ms><ns| in A
      for(uint s = 0; s < pow(2,qs_trace.size()); s++){
        uint a_out = a_m, a_in = a_n;
        for(uint q = 0; q < qs_trace.size(); q++){
          if(qbit_state(q,qs_trace.size(),s)){
            a_out += bit_int(qs_trace.at(q),qbits_old);
            a_in += bit_int(qs_trace.at(q),qbits_old);
          }
        }
        B(m,n) += A(a_out,a_in);
      }
    }
  }
  return B;
}

//--------------------------------------------------------------------------------------------
// Matrix vectors
//--------------------------------------------------------------------------------------------

mvec operator*(const MatrixXcd& G, const mvec& v){
  vector<MatrixXcd> out;
  for(uint i = 0; i < v.size(); i++){
    out.push_back(G*v.at(i));
  }
  return mvec(out);
}

