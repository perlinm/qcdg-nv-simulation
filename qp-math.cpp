#define EIGEN_USE_MKL_ALL

#include <iostream> // for standard output
#include <vector> // vector objects
#include <eigen3/Eigen/Dense> // linear algebra library
#include <eigen3/unsupported/Eigen/KroneckerProduct> // provides tensor product
#include <eigen3/unsupported/Eigen/MatrixFunctions> // provides matrix functions

#include "constants.h"
#include "qp-math.h"

using namespace std;
using namespace Eigen;

// ---------------------------------------------------------------------------------------
// Matrix functions
// ---------------------------------------------------------------------------------------

// tensor product of many matrices
MatrixXcd tp(const initializer_list<MatrixXcd>& list) {
  MatrixXcd out = I1;
  for (MatrixXcd elem: list) {
    out = tp(out,elem);
  }
  return out;
}

// remove numerical artifacts from a matrix
MatrixXcd remove_artifacts(const MatrixXcd& A, const double threshold) {
  MatrixXcd B = MatrixXcd::Zero(A.rows(),A.cols());
  for (uint m = 0; m < A.rows(); m++) {
    for (uint n = 0; n < A.cols(); n++) {
      B(m,n) += round(A(m,n).real()/threshold)*threshold;
      B(m,n) += round(A(m,n).imag()/threshold)*threshold*j;
    }
  }
  return B;
}

// get global phase of matrix
complex<double> get_phase(const MatrixXcd& A, const double threshold) {
  for (uint m = 0; m < A.rows(); m++) {
    for (uint n = 0; n < A.cols(); n++) {
      if (abs(A(m,n)) > threshold) {
        const complex<double> phase = A(m,n)/abs(A(m,n));
        return phase;
      }
    }
  }
  return 1;
}

// ---------------------------------------------------------------------------------------
// Operator manipulation
// ---------------------------------------------------------------------------------------

// generate matrix B to act A on qbits qs_act out of qbits_new
MatrixXcd act(const MatrixXcd& A, const vector<uint>& qs_act, const uint qbits_new) {
  assert(A.rows() == A.cols()); // A should be square

  if (qs_act.size() == qbits_new) {
    bool do_nothing = true;
    for (uint i = 0; i < qbits_new; i++) {
      if (i != qs_act.at(i)) {
        do_nothing = false;
        break;
      }
    }
    if (do_nothing) return A;
  }

  // number of qbits A acted on
  const uint qbits_old = qs_act.size();
  assert(qbits_old == log2(A.rows()));

  // vector of qbits we are ignoring
  vector<uint> qs_ignore;
  for (uint i = 0; i < qbits_new; i++) {
    if (!in_vector(i,qs_act)) {
      qs_ignore.push_back(i);
    }
  }

  // initialize B (output) to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries A(m,n)
  for (uint m = 0; m < A.rows(); m++) {
    for (uint n = 0; n < A.cols(); n++) {

      // get contribution of substates |m-><n-| to indices of B
      uint b_m = 0, b_n = 0;
      for (uint q = 0; q < qbits_old; q++) {
        if (qbit_state(q,qbits_old,m)) b_m += bit_int(qs_act.at(q),qbits_new);
        if (qbit_state(q,qbits_old,n)) b_n += bit_int(qs_act.at(q),qbits_new);
      }

      // loop over all elements of the form |ms><ns| in B
      for (uint s = 0; s < pow(2,qs_ignore.size()); s++) {
        uint b_out = b_m, b_in = b_n;
        for (uint q = 0; q < qs_ignore.size(); q++) {
          if (qbit_state(q,qs_ignore.size(),s)) {
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
MatrixXcd ptrace(const MatrixXcd& A, const vector<uint>& qs_trace) {
  if (qs_trace.size() == 0) return A;

  assert(A.rows() == A.cols()); // A should be square

  // number of qbits A acted on
  const uint qbits_old = log2(A.rows());
  const uint qbits_new = qbits_old - qs_trace.size();

  // vector of qbits we are keeping
  vector<uint> qs_keep;
  for (uint i = 0; i < qbits_old; i++) {
    if (!in_vector(i,qs_trace)) qs_keep.push_back(i);
  }
  assert(qbits_new == qs_keep.size());

  // initialize B (output) to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries B(m,n)
  for (uint m = 0; m < B.rows(); m++) {
    for (uint n = 0; n < B.cols(); n++) {

      // get contribution of substates |m-><n-| to indices of A
      uint a_m = 0, a_n = 0;
      for (uint q = 0; q < qbits_new; q++) {
        if (qbit_state(q,qbits_new,m)) a_m += bit_int(qs_keep.at(q),qbits_old);
        if (qbit_state(q,qbits_new,n)) a_n += bit_int(qs_keep.at(q),qbits_old);
      }

      // loop over all elements of the form |ms><ns| in A
      for (uint s = 0; s < pow(2,qs_trace.size()); s++) {
        uint a_out = a_m, a_in = a_n;
        for (uint q = 0; q < qs_trace.size(); q++) {
          if (qbit_state(q,qs_trace.size(),s)) {
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

// ---------------------------------------------------------------------------------------
// Gate decomposition and fidelity
// ---------------------------------------------------------------------------------------

// returns basis element p for an operator acting on a system with N spins
MatrixXcd U_basis_element(const uint p, const uint N) {
  const vector<MatrixXcd> spins = {I2, sx, sy, sz};
  MatrixXcd b_p = I1;
  for (uint n = 0; n < N; n++) {
    b_p = tp(b_p,spins.at(int_bit(p,2*n)+2*int_bit(p,2*n+1)));
  }
  return b_p;
}

// returns element p of a basis for operators acting on a system with N qbits
string U_basis_element_text(const uint p, const uint N) {
  const vector<char> spins = {'I','X','Y','Z'};
  stringstream stream;
  for (uint n = 0; n < N; n++) {
    stream << spins.at(int_bit(p,2*n)+2*int_bit(p,2*n+1));
  }
  return stream.str();
}

// returns matrix whose columns are basis Hamiltonians for a system of N spins
MatrixXcd U_basis_matrix(const uint N) {
  MatrixXcd out = MatrixXcd::Zero(pow(4,N),pow(4,N));
  for (uint p = 0; p < pow(4,N); p++) {
    out.col(p) = flatten(U_basis_element(p,N));
  }
  return out;
}

// decompose an operator into its basis elements
VectorXcd U_decompose(const MatrixXcd& U, const bool fast) {
  const uint N = log2(U.rows());
  if (fast) return U_basis_matrix(N).partialPivLu().solve(flatten(U));
  else return U_basis_matrix(N).fullPivHouseholderQr().solve(flatten(U));
}

// return k-th kraus operator of U when acting on a given subsystem
MatrixXcd kraus_element(const MatrixXcd& U, const vector<uint>& system, const uint k) {
  const uint total_qbits = log2(U.rows());
  const uint system_qbits = system.size();
  const uint environment_qbits = total_qbits - system_qbits;
  const uint dim_s = pow(2,system_qbits);
  const uint dim_e = pow(2,environment_qbits);
  assert(k < dim_e);

  // rearrange U so that the system qbits are separated from the environment qbits
  vector<uint> qbit_order = {};
  for (uint qq = 0; qq < total_qbits; qq++) {
    if (in_vector(qq,system)) qbit_order.push_back(qq);
  }
  for (uint qq = 0; qq < total_qbits; qq++) {
    if (!in_vector(qq,qbit_order)) qbit_order.push_back(qq);
  }

  vector<uint> qbit_permutation = {};
  for (uint qq = 0; qq < total_qbits; qq++) {
    for (uint pos = 0; pos < total_qbits; pos++) {
      if (qbit_order.at(pos) == qq) {
        qbit_permutation.push_back(pos);
      }
    }
  }
  const MatrixXcd M = act(U,qbit_permutation,total_qbits);

  MatrixXcd E_k = MatrixXcd::Zero(dim_s,dim_s);
  for (uint i = 0; i < E_k.rows(); i++) {
    for (uint j = 0; j < E_k.cols(); j++) {
      E_k(i,j) = M(i*dim_e+k,j*dim_e);
    }
  }
  return E_k;
}

// compute mean fidelity of propagator U within a given subsystem with respect to gate G
double gate_fidelity(const MatrixXcd& U, const MatrixXcd& G, const vector<uint>& system) {
  assert(U.size() >= G.size());

  const uint spins = log2(U.rows());
  vector<uint> environment = {};
  if (system.size() > 0) {
    for (uint n = 0; n < spins; n++) {
      if (!in_vector(n,system)) {
        environment.push_back(n);
      }
    }
  }

  if (environment.size() == 0) {
    // use straightforward fidelity formula
    const MatrixXcd M = G * U.adjoint();
    const uint D = M.rows();
    return real( trace(M.adjoint()*M) + trace(M)*conj(trace(M)) ) / (D*(D+1));

  } else {
    // use fidelity formula with kraus operators for action of U on the given subsystem
    const uint dim_s = pow(2,system.size());
    const uint dim_e = pow(2,environment.size());
    MatrixXcd sum_to_trace = MatrixXcd::Zero(dim_s,dim_s);
    double sum_trace_squared = 0;
    for (uint k = 0; k < dim_e; k++) {
      const MatrixXcd M_k = G.adjoint() * kraus_element(U, system, k);
      sum_to_trace += M_k.adjoint() * M_k;
      sum_trace_squared += real(trace(M_k)*conj(trace(M_k)));
    }
    return (real(trace(sum_to_trace)) + sum_trace_squared) / (dim_s * (dim_s+1));
  }
}

// ---------------------------------------------------------------------------------------
// Matrix vectors
// ---------------------------------------------------------------------------------------

mvec operator*(const MatrixXcd& G, const mvec& v) {
  vector<MatrixXcd> out;
  for (uint i = 0; i < v.size(); i++) {
    out.push_back(G*v.at(i));
  }
  return mvec(out);
}

// ---------------------------------------------------------------------------------------
// Printing methods
// ---------------------------------------------------------------------------------------

// print operator in human-readable form
void U_print(const MatrixXcd& U, const double threshold) {
  const int N = log2(U.rows());
  const VectorXcd hs = remove_artifacts(U_decompose(U), threshold);
  for (int p = 0; p < pow(4,N); p++) {
    if (abs(hs(p)) != 0) {
      cout << U_basis_element_text(p,N) << ": " << hs(p) << endl;
    }
  }
  cout << endl;
}

// print state vector in human readable form
void state_print(const MatrixXcd& psi) {
  const uint N = psi.size();
  const uint qbits = log2(N);
  for (uint n = 0; n < N; n++) {
    if (abs(psi(n)) != 0) {
      cout << "|";
      for (uint q = 0; q < qbits; q++) {
        cout << (qbit_state(q,qbits,n)?"d":"u");
      }
      cout << "> " << psi(n) << endl;
    }
  }
  cout << endl;
}

// print matrix in human readable form
void matrix_print(const MatrixXcd& M) {
  const uint qbits = log2(M.rows());
  for (uint m = 0; m < M.rows(); m++) {
    for (uint n = 0; n < M.cols(); n++) {
      if (abs(M(m,n)) != 0) {
        cout << "|";
        for (uint q = 0; q < qbits; q++) {
          cout << (qbit_state(q,qbits,m)?"d":"u");
        }
        cout << "><";
        for (uint q = 0; q < qbits; q++) {
          cout << (qbit_state(q,qbits,n)?"d":"u");
        }
        cout << "| " << M(m,n) << endl;
      }
    }
  }
  cout << endl;
}
