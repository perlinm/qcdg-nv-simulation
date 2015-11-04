#pragma once

#include <eigen3/Eigen/Dense> // linear algebra library
#include <eigen3/unsupported/Eigen/KroneckerProduct> // provides tensor product
#include <eigen3/unsupported/Eigen/MatrixFunctions> // provides matrix functions
using namespace Eigen;

// check whether value val is in vector vec
inline bool in_vector(auto val, vector<auto> vec){
  return (find(vec.begin(), vec.end(), val) != vec.end());
}

// identity matrices
const MatrixXcd I1 = MatrixXcd::Identity(1,1);
const MatrixXcd I2 = MatrixXcd::Identity(2,2);
const MatrixXcd I4 = MatrixXcd::Identity(4,4);

// return unit vector in direction of vec
inline Vector3d hat(Vector3d vec){ return vec/vec.norm(); }

// matrix functions
inline complex<double> trace(const MatrixXcd M){ return M.trace(); }
inline MatrixXcd log(const MatrixXcd M){ return M.log(); }
inline MatrixXcd exp(const MatrixXcd M){ return M.exp(); }
inline MatrixXcd sqrt(const MatrixXcd M){ return M.sqrt(); }
inline MatrixXcd pow(const MatrixXcd M, auto x){ return M.pow(x); }

// tensor product of two matrices
inline MatrixXcd tp(const MatrixXcd A, const MatrixXcd B){ return kroneckerProduct(A,B); }

// tensor product of many matrices
MatrixXcd tp(const initializer_list<MatrixXcd> list){
  MatrixXcd out = I1;
  for(MatrixXcd elem: list){
    out = tp(out,elem);
  }
  return out;
}

// remove numerical artifacts from a matrix
void remove_artifacts(MatrixXcd &A, double threshold = 1e-12){
  for(uint m = 0; m < A.rows(); m++){
    for(uint n = 0; n < A.cols(); n++){
      if(abs(A(m,n).real()) < threshold) A(m,n) -= A(m,n).real();
      if(abs(A(m,n).imag()) < threshold) A(m,n) -= A(m,n).imag()*j;
    }
  }
}

// get global phase of matrix
complex<double> get_phase(const MatrixXcd A){
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

// remove global phase from matrix
void remove_phase(MatrixXcd &A){ A *= conj(get_phase(A)); }

//--------------------------------------------------------------------------------------------
// Operator rearrangement
//--------------------------------------------------------------------------------------------

// get the n-th bit of an integer num
inline bool int_bit(uint num, uint n){
  if(pow(2,n) > num) return 0;
  else return (num >> n) & 1;
}

// get state of qbit q (of N) from enumerated state s
inline bool qbit_state(uint q, uint N, uint s){ return int_bit(s,N-1-q); }

// get integer corresponding to an 'on' state of bit p (of N)
inline uint bit_int(uint q, int N){ return pow(2,N-1-q); }

// generate matrix B to act A on qbits qs_act out of qbits_new
MatrixXcd act(const MatrixXcd A, const vector<uint> qs_act, uint qbits_new){
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
MatrixXcd ptrace(const MatrixXcd A, const vector<uint> qs_trace){
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

struct mvec{
  vector<MatrixXcd> v;

  mvec(){};
  mvec(const vector<MatrixXcd> v){ this->v = v; };
  mvec(const MatrixXcd v_mat, const Vector3d v_vec){
    for(uint i = 0; i < v_vec.size(); i++){
      v.push_back(v_vec(i)*v_mat);
    }
  };

  uint size() const { return v.size(); }
  MatrixXcd at(uint i) const { return v.at(i); }

  bool operator==(const mvec w) const {
    assert(v.size() == w.size());
    for(uint i = 0; i < v.size(); i++){
      if(v.at(i) != w.at(i)) return false;
    }
    return true;
  }
  bool operator!=(const mvec w) const { return !(*this == w); }

  mvec operator+(const mvec w) const {
    assert(v.size() == w.size());
    mvec out = v;
    for(uint i = 0; i < v.size(); i++){
      out.at(i) += w.at(i);
    }
    return out;
  }
  mvec operator-(const mvec w) const {
    assert(v.size() == w.size());
    mvec out = v;
    for(uint i = 0; i < v.size(); i++){
      out.at(i) -= w.at(i);
    }
    return out;
  }
  mvec operator*(const double s) const {
    mvec out = v;
    for(uint i = 0; i < out.size(); i++){
      out.at(i) *= s;
    }
    return out;
  }
  mvec operator/(const double s) const { return *this * (1/s); }

  MatrixXcd dot(const mvec w) const {
    assert(v.size() == w.size());
    MatrixXcd out = tp(v.at(0),w.at(0));
    for(uint i = 1; i < v.size(); i++){
      out += tp(v.at(i),w.at(i));
    }
    return out;
  }

  MatrixXcd dot(const Vector3d r) const {
    assert(v.size() == 3);
    return v.at(0)*r(0) + v.at(1)*r(1) + v.at(2)*r(2);
  }
};

inline MatrixXcd dot(mvec v, mvec w){ return v.dot(w); }
inline MatrixXcd dot(mvec v, Vector3d r){ return v.dot(r); }
inline MatrixXcd dot(Vector3d r, mvec v){ return v.dot(r); }

mvec operator*(double s, mvec v){ return v*s; }

// product between matrix and matrix vector
mvec operator*(MatrixXcd G, mvec v){
  vector<MatrixXcd> out;
  for(uint i = 0; i < v.size(); i++){
    out.push_back(G*v.at(i));
  }
  return mvec(out);
}
mvec operator*(mvec v, MatrixXcd G){
  vector<MatrixXcd> out;
  for(uint i = 0; i < v.size(); i++){
    out.push_back(v.at(i)*G);
  }
  return mvec(out);
}
