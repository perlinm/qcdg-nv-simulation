#!/usr/bin/env python3

import matrix_operations_symbolic as mo
import sympy as sym
i = complex(0,1)

########################################
# mathematical and physical constants
########################################

arccos = sym.acos

# identity matrices
def I(n): return sym.eye(n)

# pauli spin operators
sz = sym.Matrix([[1,0],[0,-1]])
sx = sym.Matrix([[0,1],[1,0]])
sy = i*sym.Matrix([[0,-1],[1,0]])

# spin-up and spin-down states
up = sym.Matrix([[1],[0]])
dn = sym.Matrix([[0],[1]])

# two-qubit states
uu = mo.tp(up,up)
ud = mo.tp(up,dn)
du = mo.tp(dn,up)
dd = mo.tp(dn,dn)

# two-qubit singlet/triplet states
S = (ud-du)/sym.sqrt(2)
T = (ud+du)/sym.sqrt(2)

# physical constants
mu0 = sym.var('\mu_0')
ge = sym.var('\gamma_e')
gc = sym.var('\gamma_C')
D = sym.var('D')
qe = sym.var('q_e')
me = sym.var('m_e')
nm = sym.var('nm')

# excited state to use for two-level NV center
ms = sym.var('m_s')

# spin matrices of NV center in restricted subspace
Sz = ms/2*(sz+I(2))
Sx = sx/sym.sqrt(2)
Sy = ms*sy/sym.sqrt(2)
