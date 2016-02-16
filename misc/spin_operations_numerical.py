#!/usr/bin/env python3

import numpy as np
import scipy.linalg as lin

import constants_numerical as const
import matrix_operations_numerical as mo
i = complex(0,1)

########################################
# gate solutions and printing
########################################

# determine Hamiltonian to generate arbitrary gates (operation time = 1)
def Hgate(G): return i*lin.logm(G)

# compute propagator for arbitrary time-independent Hamiltonians (operation time = 1)
def U(tH): return lin.expm(-i*tH)

# ordered list of basis Hamiltonians for a system of N spins
def vH(N=2,form='matrix'):

    pauli = [np.eye(2),const.sz,const.sx,const.sy]
    pauli_text = ['I','Z','X','Y']

    vec = []
    text = []
    permutation = np.zeros(N,dtype=int)
    for p in range(4**N):
        bin = '{0:b}'.format(p).zfill(2*N)
        for s in range(N):
            permutation[s] = int(bin[s*2:s*2+2],2)
        vec.append(mo.tp([ pauli[s] for s in permutation ]))
        text.append(' '.join([ pauli_text[s] for s in permutation ]))

    if form == 'matrix':
        return np.column_stack([np.ndarray.flatten(v) for v in vec])
    elif form == 'list':
        return vec
    elif form == 'text':
        return text
    else:
        return None

# decompose Hamiltonian into its basis components
def H_decompose(H):
    N = int(np.log(H.shape[0])/np.log(2))
    hs = lin.solve(vH(N),np.ndarray.flatten(H))
    return mo.remove_artifacts(hs)

# print Hamiltonian vector in human-readable form
def H_print(H,factor=1,threshold=1e-10):
    N = int(np.log(H.shape[0])/np.log(2))
    hs = mo.remove_artifacts(H_decompose(H),threshold)
    vHs = vH(N,'text')
    factor_text = 'prefactor: '
    if factor == 1:
        factor_text += 'pi'
    else:
        factor_text += 'pi/%g'%factor
    print(factor_text)
    for n in range(len(hs)):
        if abs(hs[n]) != 0:
            pretext = '%s:'%vHs[n] + (' ' if hs[n] > 0 else '')
            print(pretext,hs[n]*factor/np.pi)

# print Hamiltonian vector generating the gate G in human-readable form
def G_print(G,factor=1,threshold=1e-10): H_print(Hgate(G),factor,threshold)

# print qubit in human-readable form
def qvec_print(v):
    N = len(v)
    qbits = int(np.log(N)/np.log(2))
    for n in range(N):
        if v[n] != 0:
            s = bin(n)[2:].zfill(qbits)
            s = s.replace('0','u').replace('1','d')
            print('%s:'%s,v[n][0])

########################################
# single qubit gates
########################################

# spin propagators; Ua corresponds to a Hamiltonian # H = h s_a
def Uz(ht): return np.cos(ht)*np.eye(2)-i*np.sin(ht)*const.sz
def Ux(ht): return np.cos(ht)*np.eye(2)-i*np.sin(ht)*const.sx
def Uy(ht): return np.cos(ht)*np.eye(2)-i*np.sin(ht)*const.sy

# rotation operators
def Rz(phi): return Uz(phi/2)
def Rx(phi): return Ux(phi/2)
def Ry(phi): return Uy(phi/2)

# phase-flip, bit-flip, and Hadamard gates
Z = np.array([[1,0],[0,-1]])
X = np.array([[0,1],[1,0]])
HG = np.array([[1,1],[1,-1]])/np.sqrt(2)

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
    return np.cos(ht)*np.eye(4) - i*np.sin(ht)*mo.tp(spins)

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
    return np.array([[1,0,0,0],
                     [0,1,0,0],
                     [0,0,1,0],
                     [0,0,0,np.exp(i*phi)]])

cRdd = cR
def cRuu(phi):
    return np.array([[np.exp(i*phi),0,0,0],
                     [0,1,0,0],
                     [0,0,1,0],
                     [0,0,0,1]])
def cRud(phi):
    return np.array([[1,0,0,0],
                     [0,np.exp(i*phi),0,0],
                     [0,0,1,0],
                     [0,0,0,1]])
def cRdu(phi):
    return np.array([[1,0,0,0],
                     [0,1,0,0],
                     [0,0,np.exp(i*phi),0],
                     [0,0,0,1]])

# controlled bit-flip and controlled-NOT gates
cZ = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0,-1]])
cNOT = np.array([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0]])

# sqrt(iSWAP), iSWAP, sqrt(SWAP), and SWAP gates
riSWAP = np.array([[1, 0,            0,             0],
                   [0, 1/np.sqrt(2), i/np.sqrt(2), 0],
                   [0, i/np.sqrt(2), 1/np.sqrt(2), 0],
                   [0, 0,            0,             1]])

iSWAP = np.array([[1, 0, 0, 0],
                  [0, 0, i, 0],
                  [0, i, 0, 0],
                  [0, 0, 0, 1]])

rSWAP = np.array([[1, 0,       0,       0],
                  [0, (1+i)/2, (1-i)/2, 0],
                  [0, (1-i)/2, (1+i)/2, 0],
                  [0, 0,       0,       1]])

SWAP = np.array([[1, 0, 0, 0],
                 [0, 0, 1, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1]])

# maps (up,down) to singlet state and (down,up) to triplet-0 state
E =  np.array([[ 1, 0, 0, 1],
               [ 0, 1, 1, 0],
               [ 0,-1, 1, 0],
               [-1, 0, 0, 1]])/np.sqrt(2)
uE = lin.inv(E)

# map \ket{\u\d} -> \cos\theta\ket S + e^{i\phi}\sin\theta\ket T
def Q(theta,phi):
    alpha = np.cos(theta)
    beta = np.sin(theta)
    phase = np.exp(i*phi)
    return np.array([[ alpha+beta, 0, 0, alpha-beta],
                     [ 0, alpha+beta*phase, alpha*phase-beta, 0],
                     [ 0,-alpha+beta*phase, alpha*phase+beta, 0],
                     [-alpha+beta, 0, 0, alpha+beta]])/np.sqrt(2)

# coupled -> uncoupled basis change for two spins
uncouple = np.array([[1, 0,            0,            0],
                     [0, 1/np.sqrt(2), 1/np.sqrt(2), 0],
                     [0,-1/np.sqrt(2), 1/np.sqrt(2), 0],
                     [0, 0,            0,            1]])
couple = lin.inv(uncouple)

########################################
# ST single qubit gates
########################################

# identity on unused two-qubit subspace
Ir_ST = mo.mm(const.uu,const.uu.T) + mo.mm(const.dd,const.dd.T)

# rotation operators
def Rz_ST(phi):
    return np.exp(-i*phi/2)*mo.mm(const.S,const.S.T) \
           + np.exp(i*phi/2)*mo.mm(const.T,const.T.T) + Ir_ST

def Rx_ST(phi):
    return np.cos(phi/2)*(mo.mm(const.S,const.S.T)+mo.mm(const.T,const.T.T)) \
           - i*np.sin(phi/2)*(mo.mm(const.S,const.T.T)+mo.mm(const.T,const.S.T)) + Ir_ST

def Ry_ST(phi):
    return np.cos(phi/2)*(mo.mm(const.S,const.S.T)+mo.mm(const.T,const.T.T)) \
           + np.sin(phi/2)*(-mo.mm(const.S,const.T.T)+mo.mm(const.T,const.S.T)) + Ir_ST

# phase-flip, bit-flip, and Hadamard gates
Z_ST = mo.mm(const.S,const.S.T)-mo.mm(const.T,const.T.T)+i*Ir_ST
X_ST = mo.mm(const.S,const.T.T)+mo.mm(const.T,const.S.T)+i*Ir_ST
HG_ST = (mo.mm(const.S+const.T,const.S.T) \
         + mo.mm(const.S-const.T,const.T.T))/np.sqrt(2) + i*Ir_ST
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

