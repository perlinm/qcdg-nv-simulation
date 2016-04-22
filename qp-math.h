#pragma once

#include <eigen3/Eigen/Dense> // linear algebra library
#include <eigen3/unsupported/Eigen/KroneckerProduct> // provides tensor product
#include <eigen3/unsupported/Eigen/MatrixFunctions> // provides matrix functions
using namespace Eigen;

// floating-point modulus
inline double mod(const double number, const double modulus){
  return number - floor(number/modulus)*modulus;
}

// check whether val is in vec
inline bool in_vector(const uint val, const vector<uint>& vec){
  return (find(vec.begin(), vec.end(), val) != vec.end());
}

// dot product between vectors
inline double dot(const Vector3d& v, const Vector3d& w){ return v.dot(w); }

// return unit vector in direction of vec
inline Vector3d hat(const Vector3d& vec){ return vec.normalized(); }

// project vec onto plane orthogonal to axis
inline Vector3d project(const Vector3d& vec, const Vector3d& axis){
  return vec - vec.dot(hat(axis))*hat(axis);
}

// return a vector rotated by a given angle about a given axis
inline Vector3d rotate(const Vector3d& vec, const double angle, const Vector3d& axis){
  return (vec.dot(hat(axis))*hat(axis) +
          cos(angle) * project(vec,axis) +
          sin(angle) * hat(axis).cross(vec));
}

//--------------------------------------------------------------------------------------------
// Matrix functions
//--------------------------------------------------------------------------------------------

inline complex<double> trace(const MatrixXcd& M){ return M.trace(); }
inline MatrixXcd log(const MatrixXcd& M){ return M.log(); }
inline MatrixXcd exp(const MatrixXcd& M){ return M.exp(); }
inline MatrixXcd sqrt(const MatrixXcd& M){ return M.sqrt(); }
inline MatrixXcd pow(const MatrixXcd& M, const double x){ return M.pow(x); }

// tensor product of two matrices
inline MatrixXcd tp(const MatrixXcd& A, const MatrixXcd& B){ return kroneckerProduct(A,B); }

// tensor product of many matrices
MatrixXcd tp(const initializer_list<MatrixXcd>& list);

// remove numerical artifacts from a matrix
MatrixXcd remove_artifacts(const MatrixXcd& A, const double threshold = 1e-12);

// get global phase of matrix
complex<double> get_phase(const MatrixXcd& A, const double threshold = 1e-12);

// remove global phase from matrix
inline MatrixXcd remove_phase(const MatrixXcd& A, const double threshold = 1e-12){
  return A*conj(get_phase(A, threshold));
}

// clean up matrix for human readability
inline MatrixXcd clean(const MatrixXcd& M, double threshold = 1e-3){
  return remove_artifacts(remove_phase(M,threshold),threshold);
}

//--------------------------------------------------------------------------------------------
// Operator manipulation
//--------------------------------------------------------------------------------------------

// get the n-th bit of an integer num
inline bool int_bit(const uint num, const uint n){
  if(pow(2,n) > num) return 0;
  else return (num >> n) & 1;
}

// get state of qbit q (of N) from enumerated state s
inline bool qbit_state(const uint q, const uint N, const uint s){ return int_bit(s,N-1-q); }

// get integer corresponding to an 'on' state of bit p (of N)
inline uint bit_int(const uint q, const int N){ return pow(2,N-1-q); }

// generate matrix B to act A on qbits qs_act out of qbits_new
MatrixXcd act(const MatrixXcd& A, const vector<uint>& qs_act, const uint qbits_new);

// perform a partial trace over qbits qs_trace
MatrixXcd ptrace(const MatrixXcd& A, const vector<uint>& qs_trace);

//--------------------------------------------------------------------------------------------
// Gate decomposition and fidelity
//--------------------------------------------------------------------------------------------

// returns element p of a basis for operators acting on a system with N qubits
MatrixXcd U_basis_element(const uint p, const uint N);

// returns element p of a basis for operators acting on a system with N qubits
string U_basis_element_text(const uint p, const uint N);

// flatten matrix into a 1-D vector
inline MatrixXcd flatten(MatrixXcd M){
  M.resize(M.size(),1);
  return M;
}

// returns matrix whose columns are basis operators for a system of N spins
MatrixXcd U_basis_matrix(const uint N);

// decompose Hamiltonian into its basis elements
VectorXcd U_decompose(const MatrixXcd& U, const bool fast = false);

// compute mean fidelity of gate U with respect to G, i.e. how well U approximates G
double gate_fidelity(const MatrixXcd&  U, const MatrixXcd& G);

// compute mean fidelity of propagator acting on given nuclei
double gate_fidelity(const MatrixXcd& U, const MatrixXcd& G, const vector<uint>& nuclei);

// compute fidelity of state rho with respect to state sigma, i.e. how close rho is to sigma
inline double state_fidelity(const MatrixXcd& rho, const MatrixXcd& sigma){
  const MatrixXcd sqrt_rho = sqrt(rho);
  const double sqrt_F = abs(trace(sqrt(sqrt_rho*sigma*sqrt_rho)));
  return sqrt_F*sqrt_F;
}

//--------------------------------------------------------------------------------------------
// Matrix vectors
//--------------------------------------------------------------------------------------------

struct mvec{
  vector<MatrixXcd> v;

  mvec(){};
  mvec(const vector<MatrixXcd>& v){ this->v = v; };
  mvec(const MatrixXcd& v_mat, const Vector3d& v_vec){
    for(uint i = 0; i < v_vec.size(); i++){
      v.push_back(v_mat*v_vec(i));
    }
  }

  uint size() const { return v.size(); }
  MatrixXcd at(const uint i) const { return v.at(i); }

  // comparison operators
  bool operator==(const mvec& w) const {
    assert(v.size() == w.size());
    for(uint i = 0; i < v.size(); i++){
      if(v.at(i) != w.at(i)) return false;
    }
    return true;
  }
  bool operator!=(const mvec& w) const { return !(*this == w); }

  // addition, subtraction, and multiplication
  mvec operator+(const mvec& w) const {
    assert(v.size() == w.size());
    vector<MatrixXcd> out;
    for(uint i = 0; i < v.size(); i++){
      out.push_back(v.at(i)+w.at(i));
    }
    return out;
  }
  mvec operator-(const mvec& w) const {
    assert(v.size() == w.size());
    vector<MatrixXcd> out;
    for(uint i = 0; i < v.size(); i++){
      out.push_back(v.at(i)-w.at(i));
    }
    return out;
  }
  mvec operator*(const double s) const {
    vector<MatrixXcd> out;
    for(uint i = 0; i < v.size(); i++){
      out.push_back(v.at(i)*s);
    }
    return out;
  }
  mvec operator/(const double s) const { return *this * (1/s); }
  mvec operator/(const int s) const { return *this * (1/double(s)); }
  mvec operator*(const MatrixXcd& G) const {
    vector<MatrixXcd> out;
    for(uint i = 0; i < v.size(); i++){
      out.push_back(v.at(i)*G);
    }
    return out;
  }

  // inner product with vectors and matrix vectors
  MatrixXcd dot(const Vector3d& r) const {
    assert(v.size() == 3);
    return v.at(0)*r(0) + v.at(1)*r(1) + v.at(2)*r(2);
  }
  MatrixXcd dot(const mvec& w) const {
    assert(v.size() == w.size());
    MatrixXcd out = tp(v.at(0),w.at(0));
    for(uint i = 1; i < v.size(); i++){
      out += tp(v.at(i),w.at(i));
    }
    return out;
  }
};

inline MatrixXcd dot(const mvec& v, const mvec& w){ return v.dot(w); }
inline MatrixXcd dot(const mvec& v, const Vector3d& r){ return v.dot(r); }
inline MatrixXcd dot(const Vector3d& r, const mvec& v){ return v.dot(r); }

inline mvec operator*(const double s, mvec& v){ return v*s; }
mvec operator*(const MatrixXcd& G, const mvec& v);

//--------------------------------------------------------------------------------------------
// Common constant objects
//--------------------------------------------------------------------------------------------

// identity matrices
const Matrix<std::complex<double>,1,1> I1 = MatrixXcd::Identity(1,1);
const Matrix2cd I2 = Matrix2cd::Identity();
const Matrix4cd I4 = Matrix4cd::Identity();

// spin up/down state vectors
const Vector2cd up = (Vector2cd() << 1,0).finished();
const Vector2cd dn = (Vector2cd() << 0,1).finished();

// two qbit basis states
const Vector4cd uu = tp(up,up);
const Vector4cd ud = tp(up,dn);
const VectorXcd du = tp(dn,up);
const Vector4cd dd = tp(dn,dn);

// pauli spin matrices
const Matrix2cd sx = up*dn.adjoint() + dn*up.adjoint();
const Matrix2cd sy = j*(-up*dn.adjoint() + dn*up.adjoint());
const Matrix2cd sz = up*up.adjoint() - dn*dn.adjoint();

//--------------------------------------------------------------------------------------------
// Printing methods
//--------------------------------------------------------------------------------------------

// print operator in human-readable form
void U_print(const MatrixXcd& U, const double threshold = 1e-3);

// print state vector in human readable form
void state_print(const MatrixXcd& psi);

// print matrix in human readable form
void matrix_print(const MatrixXcd& M);
