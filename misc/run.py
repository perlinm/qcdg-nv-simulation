#!/usr/bin/env python3

import sys
from scipy.linalg import expm,logm,sqrtm

if 'sym' in sys.argv:
    import pylab as pl
    from sympy import *
    from constants_symbolic import *
    from matrix_operations_symbolic import *
    from spin_operations_symbolic import *
    matrix = Matrix
else:
    from pylab import *
    from constants_numerical import *
    from matrix_operations_numerical import *
    from spin_operations_numerical import *
    import numpy as np
    matrix = array

    np.set_printoptions(linewidth = 95)




phi_over_2pi = float(sys.argv[1])
U = expm(-i*phi_over_2pi*2*pi*tp(sz,sx))
print(clean(U))


exit()

coupled = True

def pp(psi):
    psi_NVST = mm(act(couple,3,1,2),psi) if coupled else psi
    qvec_print(remove_artifacts(psi_NVST))
    print()

for psi in [tp(up,S),tp(up,T),tp(dn,S),tp(dn,T)]:
    pp(psi)
    pp(mm(SWAP_NVST,psi))
    print()
