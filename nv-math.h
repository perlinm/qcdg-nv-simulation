#pragma once

#include <eigen3/Eigen/Dense> // linear algebra library
#include <eigen3/unsupported/Eigen/KroneckerProduct> // provides tensor product
#include <eigen3/unsupported/Eigen/MatrixFunctions> // provides matrix functions
using namespace Eigen;

const complex<double> j(0.,1.);
const double pi = M_PI;

//--------------------------------------------------------------------------------------------
// Physical constants
//--------------------------------------------------------------------------------------------

// physical constants in SI
const double mu0_SI = 4*pi*1e-7; // magnetic constant (tesla meters / amp)
const double c_SI = 299792458; // speed of light (meters / second)
const double hbar_SI = 6.582119514e-16; // reduced Planck constant (eV / second)
const double qe_SI = 1.60217646e-19; // charge of electron (coulombs)
const double ge_SI = 1.760859708e11; // gyromagnetic ratio of NV electron (Hz/tesla)
const double gC13_SI = 67.28284e6; // gyromagnetic ratio of C-13 (Hz/tesla)
const double gN15_SI = -27.116e6; // gyromagnetic ratio of N-15 (Hz/tesla)
const double gH1_SI = 267.513e6; // gyromagnetic ratio of H-1 (Hz/tesla)

// physical constants in natural units: e_0 = mu_0 = c = hbar = 1; basic units are s and Hz
const double alpha = 1/137.035999074; // fine structure constant
const double qe = sqrt(4*pi*alpha);; // unit electric charge
const double me = 510998.928/hbar_SI; // mass of electron (Hz)
const double ge = qe/(2*me) * 2.0023193043617; // gyromagnetic ratio of electron
const double gC13 = gC13_SI * ge/ge_SI; // gyromagnetic ratio of C-13
const double gN15 = gN15_SI * ge/ge_SI; // gyromagnetic ratio of N-15
const double gH1 = gH1_SI * ge/ge_SI; // gyromagnetic ratio of H-1

// SI values in natural units
const double meter = 1 / c_SI; // one meter (s)
const double nm = 1e-9*meter; // one nanometer (s)
const double coulomb = qe / qe_SI; // one coulomb
const double tesla = ge_SI / ge; // one tesla (Hz)
const double volt = 1/(qe*hbar_SI); // one volt (Hz)

const double D = 2*pi*2.87e9; // NV center zero field splitting energy (Hz)
const double a0 = 0.35668 * nm; // diamond lattice parameter at 300 K
const double c13_natural_abundance = 0.0107; // natural isotopic abundance of C-13

// diamond lattice vectors scaled to a0
const Vector3d ao = (Vector3d() << 1,1,1).finished()/4;
const Vector3d a1 = (Vector3d() << 0,1,1).finished()/2;
const Vector3d a2 = (Vector3d() << 1,0,1).finished()/2;
const Vector3d a3 = (Vector3d() << 1,1,0).finished()/2;

// vector of lattice sites in a diamond unit cell
const vector<Vector3d> cell_sites { Vector3d::Zero(), a1, a2, a3, ao, ao+a1, ao+a2, ao+a3 };

//--------------------------------------------------------------------------------------------
// Useful mathematical constants and methods
//--------------------------------------------------------------------------------------------

// identity matrices
const MatrixXcd I1 = Matrix<complex<double>,1,1>::Identity();
const MatrixXcd I2 = Matrix<complex<double>,2,2>::Identity();
const MatrixXcd I4 = Matrix<complex<double>,4,4>::Identity();

// return unit vector in direction of vec
inline Vector3d normed(Vector3d vec){ return vec/vec.norm(); }

// matrix functions
inline MatrixXcd log(MatrixXcd M){ return M.log(); }
inline MatrixXcd exp(MatrixXcd M){ return M.exp(); }
inline MatrixXcd sqrt(MatrixXcd M){ return M.sqrt(); }
inline MatrixXcd pow(MatrixXcd M, auto x){ return M.pow(x); }

// tensor product of matrices in list
MatrixXcd tp(initializer_list<MatrixXcd> list){
  MatrixXcd out = (Matrix<complex<double>,1,1>() << 1).finished();
  for(MatrixXcd elem: list){
    out = kroneckerProduct(out,elem).eval();
  }
  return out;
}

// inner product of two spin vectors
inline MatrixXcd operator*(const vector<MatrixXcd> v, const vector<MatrixXcd> w){
  return tp({v.at(0),v.at(0)}) + tp({v.at(1),v.at(1)}) + tp({v.at(2),v.at(2)});
}

// inner product of a spin vector and a spacial vector
inline MatrixXcd operator*(const vector<MatrixXcd> v, const Vector3d r){
  return v.at(0)*r(0) + v.at(1)*r(1) + v.at(2)*r(2);
}
inline MatrixXcd operator*(const Vector3d r, const vector<MatrixXcd> v){ return v*r; }

// remove numerical artifacts from a matrix
void remove_artifacts(MatrixXcd &A, double threshold = 1e-12){
  for(int m = 0; m < A.rows(); m++){
    for(int n = 0; n < A.cols(); n++){
      if(abs(A(m,n).real()) < threshold) A(m,n) -= A(m,n).real();
      if(abs(A(m,n).imag()) < threshold) A(m,n) -= A(m,n).imag()*j;
    }
  }
}

// get global phase of matrix
complex<double> get_phase(const MatrixXcd A){
  for(int m = 0; m < A.rows(); m++){
    for(int n = 0; n < A.cols(); n++){
      if(abs(A(m,n)) != 0){
        complex<double> phase = A(m,n)/abs(A(m,n));
        if(phase.real() < 0) phase = -phase;
        return phase;
      }
    }
  }
  return 1.+j-j;
}

// remove global phase from matrix
void remove_phase(MatrixXcd &A){
  A *= conj(get_phase(A));
}


//--------------------------------------------------------------------------------------------
// Spin matrices and vectors
//--------------------------------------------------------------------------------------------

// pauli spin matrices
const MatrixXcd st = I2;
const MatrixXcd sx = (Matrix2cd() << 0, 1, 1, 0).finished();
const MatrixXcd sy = (Matrix2cd() << 0, -j, j, 0).finished();
const MatrixXcd sz = (Matrix2cd() << 1, 0, 0, -1).finished();

// spin-1 matrices of the NV center in the {|0>,|ms>} subspace
inline MatrixXcd Sx(int ms) { return sx/sqrt(2); }
inline MatrixXcd Sy(int ms) { return ms*sy/sqrt(2); }
inline MatrixXcd Sz(int ms) { return ms/2*(sz-st); }

// C-13 nuclear spin vector
const vector<MatrixXcd> S_C13 = { sx/2, sy/2, sz/2 };

// spin up/down state vectors
const MatrixXcd up = (Vector2cd() << 1, 0).finished();
const MatrixXcd dn = (Vector2cd() << 0, 1).finished();

// two qbit basis states
const MatrixXcd uu = tp({up,up});
const MatrixXcd ud = tp({up,dn});
const MatrixXcd du = tp({dn,up});
const MatrixXcd dd = tp({dn,dn});

// two qbit singlet-triplet states
const MatrixXcd S = (ud-du)/sqrt(2);
const MatrixXcd T = (ud+du)/sqrt(2);


//--------------------------------------------------------------------------------------------
// Qbit and logic gate methods
//--------------------------------------------------------------------------------------------

// determine Hamiltonian to generate gate G (with |H|*t = 1)
inline MatrixXcd Hgate(const MatrixXcd G){ return j*log(G); }

// compute propagator for time-independent Hamiltonian H (with |H|*t = 1)
inline MatrixXcd U(const MatrixXcd H){ return exp(-j*H); }

// get the n-th bit of an integer num
inline bool int_bit(int num, int n){
  if(pow(2,n) > num) return 0;
  else return (num >> n) & 1;
}

// get state of qbit p (of N) from enumerated state s
inline bool qbit_state(int p, int N, int s){ return int_bit(s,N-1-p); }

// generate matrix B to act A on qbits qs
MatrixXcd act(const MatrixXcd A, int qbits_new, initializer_list<int> qs_list){

  // number of qbits to act on
  const int qbits_old = qs_list.size();
  assert(qbits_old == log2(A.cols()));

  // make vector of qs
  vector<int> qs(qbits_old);
  int q = 0;
  for(int elem: qs_list){
    qs[q] = elem;
    q++;
  }

  // initialize B to the zero matrix
  MatrixXcd B = MatrixXcd::Zero(pow(2,qbits_new),pow(2,qbits_new));

  // loop over all entries B(m,n)
  for(int m = 0; m < B.rows(); m++){
    for(int n = 0; n < B.cols(); n++){

      // initially assume we are setting B(m,n) to an entry of A
      bool set_this_entry = true;

      // check input/output states of unused qbits
      if(qbits_new > qbits_old){
        for(int p = 0; p < qbits_new; p++){
          bool p_in_qs = false;
          for(int q = 0; q < qbits_old; q++){
            if(p == qs[q]){
              p_in_qs = true;
              break;
            }
          }
          if(!p_in_qs){
            // only set B(m,n) if input/output states of unused qbits are equal
            set_this_entry = (qbit_state(p,qbits_new,m) == qbit_state(p,qbits_new,n));
            if(!set_this_entry) break;
          }
        }
      }

      // if input/output states of unused qbits are equal
      if(set_this_entry){
        // determine which entry of A to use
        int a_out = 0, a_in = 0;
        for(int q = 0; q < qbits_old; q++){
          if(qbit_state(qs[q],qbits_new,m)) a_out += pow(2,qbits_old-1-q);
          if(qbit_state(qs[q],qbits_new,n)) a_in += pow(2,qbits_old-1-q);
        }
        B(m,n) = A(a_out,a_in);
      }
    }
  }
  return B;
}

/**************************************/
// single qbit gates
/**************************************/

// spin propagators; Ua corresponds to a Hamiltonian # H = h s_a
inline MatrixXcd Ux(double ht){ return cos(ht)*I2 - j*sin(ht)*sx; }
inline MatrixXcd Uy(double ht){ return cos(ht)*I2 - j*sin(ht)*sy; }
inline MatrixXcd Uz(double ht){ return cos(ht)*I2 - j*sin(ht)*sz; }

// rotation operators
inline MatrixXcd Rx(double phi){ return Ux(phi/2); }
inline MatrixXcd Ry(double phi){ return Uy(phi/2); }
inline MatrixXcd Rz(double phi){ return Uz(phi/2); }

// phase-flip, bit-flip, and Hadamard gates
const MatrixXcd Z = sz;
const MatrixXcd X = sx;
const MatrixXcd HG = (X+Z)/sqrt(2);

/**************************************/
// two qbit gates
/**************************************/

// spin coupling propagators; Uab corresponds to a Hamiltonian H = h s_a^0 s_b^1
inline MatrixXcd Uxx(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sx,sx}); }
inline MatrixXcd Uxy(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sx,sy}); }
inline MatrixXcd Uxz(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sx,sz}); }
inline MatrixXcd Uyx(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sy,sx}); }
inline MatrixXcd Uyy(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sy,sy}); }
inline MatrixXcd Uyz(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sy,sz}); }
inline MatrixXcd Uzx(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sz,sx}); }
inline MatrixXcd Uzy(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sz,sy}); }
inline MatrixXcd Uzz(double ht){ return cos(ht)*I4 - j*sin(ht)*tp({sz,sz}); }

// controlled phase rotations
inline MatrixXcd cRuu(double phi){
  return (Matrix4cd() << exp(j*phi),0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1).finished();
}
inline MatrixXcd cRud(double phi){
  return (Matrix4cd() << 1,0,0,0, 0,exp(j*phi),0,0, 0,0,1,0, 0,0,0,1).finished();
}
inline MatrixXcd cRdu(double phi){
  return (Matrix4cd() << 1,0,0,0, 0,1,0,0, 0,0,exp(j*phi),0, 0,0,0,1).finished();
}
inline MatrixXcd cRdd(double phi){
  return (Matrix4cd() << 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,exp(j*phi)).finished();
}
inline MatrixXcd cR(double phi){ return cRdd(phi); }

// controlled bit-flip and controlled-NOT gates
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

/**************************************/
// singlet-triplet qbit gates
/**************************************/

// identity on unused two-qbit subspace
const MatrixXcd unused_ST_I = uu*uu.transpose() + dd*dd.transpose();

// rotation operators
inline MatrixXcd Rx_ST(double phi){
  return (cos(phi/2.)*(S*S.transpose()+T*T.transpose())
          - j*sin(phi/2.)*(S*T.transpose()+T*S.transpose()) + unused_ST_I);
}
inline MatrixXcd Ry_ST(double phi){
  return (cos(phi/2)*(S*S.transpose()+T*T.transpose())
          + sin(phi/2)*(-S*T.transpose()+T*S.transpose()) + unused_ST_I);
}
inline MatrixXcd Rz_ST(double phi){
  return (exp(-j*phi/2.)*S*S.transpose() +
          exp( j*phi/2.)*T*T.transpose() + unused_ST_I);
}

// phase-flip, bit-flip, and Hadamard gates
const MatrixXcd Z_ST = S*S.transpose() - T*T.transpose() + j*unused_ST_I;
const MatrixXcd X_ST = S*T.transpose() + T*S.transpose() + j*unused_ST_I;
const MatrixXcd HG_ST = ((S+T)*S.transpose()+(S-T)*T.transpose())/sqrt(2.) + j*unused_ST_I;

/**************************************/
// NV / ST two qbit gates
/**************************************/

// operation to store ST state into NV spin
const MatrixXcd R_NVST = act(X*HG,3,{0})*act(SWAP,3,{0,1});

// SWAP between NV and ST states
const MatrixXcd SWAP_NVST =
  act(E,3,{1,2}) * act(cNOT,3,{1,2}) * act(X*HG,3,{0}) *
  act(cNOT,3,{0,2}) * act(SWAP,3,{0,1});


//--------------------------------------------------------------------------------------------
// Printing and debugging methods
//--------------------------------------------------------------------------------------------

// clean up matrix for human readability
inline void clean(MatrixXcd &M){
  remove_artifacts(M);
  remove_phase(M);
  remove_artifacts(M);
}

// returns basis element p for Hamiltonians of a system with N spins
MatrixXcd H_basis_element(int p, int N){
  MatrixXcd spins[4] = {st,sx,sy,sz};
  MatrixXcd b_p = spins[int_bit(p,0)+2*int_bit(p,1)];
  for(int n = 1; n < N; n++){
    b_p = tp({b_p,spins[int_bit(p,2*n)+2*int_bit(p,2*n+1)]});
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
  MatrixXcd spins[4] = {st,sx,sy,sz};
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
        cout << (qbit_state(n,qbits,q)?"d":"u");
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


//--------------------------------------------------------------------------------------------
// Constants and methods specific to the simulated system
//--------------------------------------------------------------------------------------------

// positions of nitrogen atom and vacancy
const Vector3d n_pos = Vector3d::Zero();
const Vector3d e_pos = ao;

// unit vectors along lattice directions
const Vector3d zhat = (Vector3d() << 1,1,1).finished()/sqrt(3); // direction from N to V site
const Vector3d xhat = (Vector3d() << 2,-1,-1).finished()/sqrt(6);
const Vector3d yhat = (Vector3d() << 0,1,-1).finished()/sqrt(2);

// hyperfine field at pos with C-13 atom; in the frame of H_NV^GS
inline Vector3d A(const Vector3d c13_pos){
  Vector3d r = c13_pos - e_pos;
  return ge*gC13/(4*pi*pow(a0*r.norm(),3)) * (zhat - 3*(zhat.dot(r))*r);
}

// spin-spin Hamiltonian between vacancy and C-13 atom; in the frame of H_NV^GS
inline MatrixXcd H_hf_j(int ms, const Vector3d c13_pos){
  return Sz(ms)*(A(c13_pos)*S_C13);
}

// nuclear Zeeman Hamiltonian for C-13 atoms
inline MatrixXcd H_nZ_j(const Vector3d B){ return -gC13*B*S_C13; }

// spin-spin Hamiltonian between C-13 atoms
inline MatrixXcd H_nn_jk(const Vector3d pos_j, const Vector3d pos_k){
  Vector3d r = pos_k - pos_j;
  return gC13*gC13/(4*pi*pow(a0*r.norm(),3)) * (S_C13*S_C13 - 3*(S_C13*r)*(S_C13*r));
}

