#!/usr/bin/env python3

import numpy as np
import matrix_operations_numerical as mo
i = complex(0,1)

# identity matrices
def I(n): return np.eye(n)

# pauli spin operators
st = I(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = i*np.array([[0,-1],[1,0]])

# spin-up and spin-down states
up = np.vstack([1,0])
dn = np.vstack([0,1])

# two-qubit basis states
uu = mo.tp(up,up)
ud = mo.tp(up,dn)
du = mo.tp(dn,up)
dd = mo.tp(dn,dn)

# two-qubit singlet/triplet states
S = (ud-du)/np.sqrt(2)
T = (ud+du)/np.sqrt(2)

# physical constants in SI
mu0_SI = 4*np.pi*1e-7 # magnetic constant (tesla meters / amp)
c_SI = 299792458 # speed of light (meters / second)
hbar_SI = 6.582119514e-16 # reduced Planck constant (eV / second)
qe_SI = 1.60217646e-19 # charge of electron (coulombs)
ge_SI = 1.760859708e11 # gyromagnetic ratio of NV electron (Hz/tesla)
gc_SI = 67.28284e6 # gyromagnetic ratio of C-13 (Hz/tesla)

# physical constants in natural units (e0 = mu0 = c = hbar = 1); all values in seconds or Hz

D = 2*np.pi*2.87e9 # NV center zero field splitting energy (Hz)

alpha = 1/137.035999074 # fine structure constant
qe = np.sqrt(4*np.pi*alpha) # unit electric charge
me = 510998.928/hbar_SI # mass of electron (Hz)
ge = qe/(2*me) * 2.0023193043617 # gyromagnetic ratio of electron
gc = gc_SI * ge/ge_SI # gyromagnetic ratio of C-13

# SI values in natural units
kHz = 1000 # (Hz)
meter = 1 / c_SI # one meter (s)
nm = 1e-9*meter # one nanometer (s)
coulomb = qe / qe_SI # one coulomb
tesla = ge_SI / ge # one tesla (Hz)
gauss = tesla*1e-4 # one gauss (Hz)
volt = 1/(qe*hbar_SI) # one volt (Hz)

nd = 0.15445 * nm # nearest neighbor distance in diamond at a temperature of 300 K
a0 = 0.35668 * nm # lattice parameter (unit cell side length) at a temperature of 300 K

ms = 1 # excited state to use for two-level NV center; can be +/- 1

# spin matrices of NV center in restricted subspace
Sz = ms/2*(sz+I(2))
Sx = sx/np.sqrt(2)
Sy = ms*sy/np.sqrt(2)

