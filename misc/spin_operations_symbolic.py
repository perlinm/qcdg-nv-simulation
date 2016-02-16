#!/usr/bin/env python3

import sympy as sym
import scipy.linalg as lin

import constants_symbolic as const
import matrix_operations_symbolic as mo
i = complex(0,1)

########################################
# printing
########################################

# print qubit in human-readable form
def qvec_print(v):
    N = len(v)
    qbits = int(sym.log(N)/sym.log(2))
    for n in range(N):
        if v[n] != 0:
            s = bin(n)[2:].zfill(qbits)
            s = s.replace('0','u').replace('1','d')
            print('%s:'%s,v[n])

########################################
# single qubit gates
########################################

# spin propagators; Ua corresponds to a Hamiltonian # H = h s_a
def Uz(ht): return sym.cos(ht)*sym.eye(2)-i*sym.sin(ht)*const.sz
def Ux(ht): return sym.cos(ht)*sym.eye(2)-i*sym.sin(ht)*const.sx
def Uy(ht): return sym.cos(ht)*sym.eye(2)-i*sym.sin(ht)*const.sy

# rotation operators
def Rz(phi): return Uz(phi/2)
def Rx(phi): return Ux(phi/2)
def Ry(phi): return Uy(phi/2)

# phase-flip, bit-flip, and Hadamard gates
Z = sym.Matrix([[1,0],[0,-1]])
X = sym.Matrix([[0,1],[1,0]])
HG = sym.Matrix([[1,1],[1,-1]])/sym.sqrt(2)

########################################
# two qubit gates
########################################

# propagators coupling two spins; Uab corresponds to a Hamiltonian H = h s_a^1 s_b^2
def U_spins(ht,axes):
    spins = []
    for n in range(2):
        if axes[n] == 'z':
            spin = const.sz
        elif axes[n] == 'x':
            spin = const.sx
        elif axes[n] == 'y':
            spin = const.sy
        spins.append(spin)
    return sym.cos(ht)*sym.eye(4) - i*sym.sin(ht)*mo.tp(spins)

def Uzz(ht): return U_spins(ht,'zz')
def Uxx(ht): return U_spins(ht,'xx')
def Uyy(ht): return U_spins(ht,'yy')
def Uxy(ht): return U_spins(ht,'xy')
def Uyx(ht): return U_spins(ht,'yx')
def Uyz(ht): return U_spins(ht,'yz')
def Uzy(ht): return U_spins(ht,'zy')
def Uzx(ht): return U_spins(ht,'zx')
def Uxz(ht): return U_spins(ht,'xz')

# controlled phase rotations
def cR(phi):
    R = mo.tp(Rz(phi/2),Rz(phi/2))
    return mo.mm(R,Uzz(-phi/4))

cRdd = cR
def cRuu(phi):
    R = mo.tp(Rz(-phi/2),Rz(-phi/2))
    out = mo.mm(R,Uzz(-phi/4))
    return out
def cRud(phi):
    R = mo.tp(Rz(-phi/2),Rz(phi/2))
    return mo.mm(R,Uzz(phi/4))
def cRdu(phi):
    R = mo.tp(Rz(phi/2),Rz(-phi/2))
    return mo.mm(R,Uzz(phi/4))

# controlled bit-flip and controlled-NOT gates
cZ = sym.Matrix([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0],
                 [0, 0, 0,-1]])
cNOT = sym.Matrix([[1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 1],
                   [0, 0, 1, 0]])

# sqrt(iSWAP), iSWAP, sqrt(SWAP), and SWAP gates
riSWAP = sym.Matrix([[1, 0,             0,             0],
                     [0, 1/sym.sqrt(2), i/sym.sqrt(2), 0],
                     [0, i/sym.sqrt(2), 1/sym.sqrt(2), 0],
                     [0, 0,             0,             1]])

iSWAP = sym.Matrix([[1, 0, 0, 0],
                    [0, 0, i, 0],
                    [0, i, 0, 0],
                    [0, 0, 0, 1]])

rSWAP = sym.Matrix([[1, 0,       0,       0],
                    [0, (1+i)/2, (1-i)/2, 0],
                    [0, (1-i)/2, (1+i)/2, 0],
                    [0, 0,       0,       1]])

SWAP = sym.Matrix([[1, 0, 0, 0],
                   [0, 0, 1, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 1]])

# maps (up,down) to singlet state and (down,up) to triplet-0 state
E =  sym.Matrix([[ 1, 0, 0, 1],
                 [ 0, 1, 1, 0],
                 [ 0,-1, 1, 0],
                 [-1, 0, 0, 1]])/sym.sqrt(2)
uE = E.inv()

# map \ket{\u\d} -> \cos\theta\ket S + e^{i\phi}\sin\theta\ket T
def Q(theta,phi):
    alpha = sym.cos(theta)
    beta = sym.sin(theta)
    phase = sym.exp(i*phi)
    return sym.Matrix([[ alpha+beta, 0, 0, alpha-beta],
                       [ 0, alpha+beta*phase, alpha*phase-beta, 0],
                       [ 0,-alpha+beta*phase, alpha*phase+beta, 0],
                       [-alpha+beta, 0, 0, alpha+beta]])/sym.sqrt(2)

# coupled -> uncoupled basis change for two spins
uncouple = sym.Matrix([[1, 0,             0,             0],
                       [0, 1/sym.sqrt(2), 1/sym.sqrt(2), 0],
                       [0,-1/sym.sqrt(2), 1/sym.sqrt(2), 0],
                       [0, 0,             0,             1]])
couple = uncouple.inv()

########################################
# ST single qubit gates
########################################

# identity on unused two-qubit subspace
Ir_ST = mo.mm(const.uu,const.uu.T) + mo.mm(const.dd,const.dd.T)


# rotation operators
def Rz_ST(phi):
    return sym.exp(-i*phi/2)*mo.mm(const.S,const.S.T) \
           + sym.exp(i*phi/2)*mo.mm(const.T,const.T.T) + Ir_ST

def Rx_ST(phi):
    return sym.cos(phi/2)*(mo.mm(const.S,const.S.T)+mo.mm(const.T,const.T.T)) \
           - i*sym.sin(phi/2)*(mo.mm(const.S,const.T.T)+mo.mm(const.T,const.S.T)) + Ir_ST

def Ry_ST(phi):
    return sym.cos(phi/2)*(mo.mm(const.S,const.S.T)+mo.mm(const.T,const.T.T)) \
           + sym.sin(phi/2)*(-mo.mm(const.S,const.T.T)+mo.mm(const.T,const.S.T)) + Ir_ST

# phase-flip, bit-flip, and Hadamard gates
Z_ST = mo.mm(const.S,const.S.T)-mo.mm(const.T,const.T.T)+i*Ir_ST
X_ST = mo.mm(const.S,const.T.T)+mo.mm(const.T,const.S.T)+i*Ir_ST
HG_ST = (mo.mm(const.S+const.T,const.S.T) \
         + mo.mm(const.S-const.T,const.T.T))/sym.sqrt(2) + i*Ir_ST
# Z_ST = SWAP
# X_ST = mo.act(Z,2,0)
# HG_ST = mo.mm(mo.act(Z,2,0),uE)

########################################
# NV / ST two qubit gates
########################################

R_NVST = mo.mm(mo.act(mo.mm(X,HG),3,0),mo.act(SWAP,3,0,1))

cNOT_NVST = mo.act(cZ,3,0,1)
cNOT_STNV = mo.mm(mo.act(E,3,1,2),mo.act(cNOT,3,1,0),mo.act(uE,3,1,2))
SWAP_NVST = mo.mm(cNOT_NVST,cNOT_STNV,cNOT_NVST)
