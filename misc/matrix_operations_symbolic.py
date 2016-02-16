#!/usr/bin/env python3

import numpy as np
import sympy as sym

i = complex(0,1)

# hermitian conjugate
def herm(M): return M.H

# matrix inverse
def inv(M): return M.inv()

# indices of matrix
def indices(M):
    return [ (m,n) for n in range(M.shape[1]) for m in range(M.shape[0]) ]

# commutator [X,Y]
def comm(X,Y): return X*Y-Y*X

# tensor product
def tp(*args):
    if len(args) == 1 and type(args[0]) is list:
        args = args[0]
    out = np.array(args[0])
    for i in range(len(args)-1):
        out = np.kron(out,np.array(args[i+1]))
    return sym.Matrix(out)

# matrix multiplication
def mm(*args):
    if len(args) == 1 and type(args[0]) is list:
        args = args[0]
    out = sym.simplify(args[-1])
    for i in reversed(range(len(args)-1)):
        out = sym.simplify(args[i]*out)
    return out

# matrix multiplication with object evaluation
def emm(*args):
    if len(args) == 1 and type(args[0]) is list:
        args = args[0]
    args = [ arg.evalf() for arg in args ]
    out = args[-1]
    for i in reversed(range(len(args)-1)):
        out = (args[i]*out).evalf()
    return out

# get/remove global phase from matrix
def get_phase(M):
    phase = 1
    for ind in indices(M):
        val = sym.simplify(M[ind])
        if sym.Abs(val) != 0:
            phase = val/sym.Abs(sym.expand(val))
            phase = phase.expand(complex=True)
            break
    return sym.simplify(phase)

def remove_phase(M):
    return sym.simplify(sym.expand(M*sym.conjugate(get_phase(M))))

# clean up matrix for printing
def clean(M):
    out = M.evalf()
    out = remove_phase(out)
    return out

# generate matrix B to act A on qubits ns (out of N qubits total)
def act(A,N,*ns):
    if len(ns) == 1 and type(ns[0]) is list:
        ns = ns[0]
    B = sym.zeros(2**N)
    for ind in indices(B):
        s_out = bin(ind[0])[2:].zfill(N)
        rest_out = [ s_out[n] for n in range(N) if n not in ns ]

        s_in = bin(ind[1])[2:].zfill(N)
        rest_in = [ s_in[n] for n in range(N) if n not in ns ]

        if rest_out == rest_in:
            ns_out = ''.join([ s_out[n] for n in ns ])
            ns_in = ''.join([ s_in[n] for n in ns ])
            B[ind] = A[int(ns_out,2),int(ns_in,2)]

    return B

# convert object to string for printing in latex
def to_latex(obj,keep_phase=True):
    # remove phase and format it for printing
    phase_string = ''
    try:
        if not keep_phase:
            phase = get_phase(obj)
            obj = remove_phase(obj)
            if phase != 1:
                phase_string = to_latex(phase)+'\n'
    except:
        None
    # convert object to latex text and reformat
    string = sym.latex(obj)
    # insert my latex aliases
    string = string.replace(r'\left[\begin{matrix}',r'\m{')
    string = string.replace(r'\end{matrix}\right]','}')
    string = string.replace(r'\left(','\p{')
    string = string.replace(r'\left (','\p{')
    string = string.replace(r'\right)','}')
    string = string.replace(r'\right )','}')
    string = string.replace(r'\frac',r'\f')
    # deal with stray and unnecessary '1.0's
    if string[:4] in ['1.0 ','1.0*'] :
        string = string[4:]
    string = string.replace('1.0 &','1 &')
    string = string.replace(' 1.0 ',' ')
    string = string.replace('{1.0 ','{ ')
    string = string.replace('1.0}','1}')
    string = string.replace(r'1.0\\',r'1\\')
    string = string.replace(r'\\1.0 ',r'\\ ')
    # convert some common recognized numers
    string = string.replace(r'0.5 \sqrt{2}',r'\f1{\sqrt2}')
    string = string.replace(r'\sqrt{2} \p{0.5 - 0.5 i}',r'e^{i\pi/4}')
    string = string.replace(r'\f1{\sqrt2} - \f1{\sqrt2} i',r'e^{-i\pi/4}')
    # miscellaneous
    string = string.replace(r'\\',' \\\\\n ')
    string = string.replace('- \\','-\\')
    return phase_string+string

