#!/usr/bin/env python3

import numpy as np

i = complex(0,1)

# hermitian conjugate
def herm(M):
    return np.array(np.matrix(M).H)

# indices of matrix
def indices(M):
    try:
        return [ (m,n) for n in range(M.shape[1]) for m in range(M.shape[0]) ]
    except:
        return range(M.shape[0])

# commutator [X,Y]
def comm(X,Y): return np.dot(X,Y)-np.dot(Y,X)

# tensor product
def tp(*args):
    if len(args) == 1 and type(args[0]) is list:
        args = args[0]
    out = args[0]
    for i in range(len(args)-1):
        out = np.kron(out,args[i+1])
    return out

# matrix multiplication
def mm(*args):
    if len(args) == 1 and type(args[0]) is list:
        args = args[0]
    out = args[-1]
    for i in reversed(range(len(args)-1)):
        out = np.dot(args[i],out)
    return out

# remove numerical artifacts from matrix
def remove_artifacts(M, threshold=1e-12):
    for ind in indices(M):
        if abs(np.real(M[ind])) < threshold:
            M[ind] -= np.real(M[ind])
        if abs(np.imag(M[ind])) < threshold:
            M[ind] = np.real(M[ind])
    if abs(np.imag(M)).max() == 0:
        M = np.real(M)
    return M

# get/remove global phase from matrix
def get_phase(M):
    phase = 1
    for ind in indices(M):
        if abs(M[ind]) != 0:
            phase = M[ind]/abs(M[ind])
            break
    return phase

def remove_phase(M):
    return M*np.conj(get_phase(M))

# clean up matrix for printing
def clean(M, threshold=1e-12):
    out = remove_artifacts(M,threshold)
    out = remove_phase(out)
    out = remove_artifacts(out,threshold)
    if abs(np.imag(out)).max() == 0:
        out = np.real(out)
    return out

# generate matrix B to act A on qubits ns (out of N qubits total)
def act(A,N,*ns):
    if len(ns) == 1 and type(ns[0]) is list:
        ns = ns[0]
    B = i*np.zeros((2**N,2**N))
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
