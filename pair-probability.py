#!/usr/bin/env python3

import sys, math, itertools, functools
import numpy as np

if len(sys.argv) not in [2,3]:
    print("usage: {} hyperfine_cutoff_in_kHz [c13_abundance]".format(sys.argv[0]))
    exit(1)

hyperfine_cutoff = float(sys.argv[1])*1000  # cutoff for hyperfine field strength; Hz
try:
    c13_abundance = float(sys.argv[2])
except:
    c13_abundance = 0.0107

# return the sum/product of all elements in a list
def sum(list):
    return functools.reduce(lambda x,y: x+y, list)

def product(list):
    return functools.reduce(lambda x,y: x*y, list)

# physical constants in SI
c_SI = 299792458 # speed of light (meters/second)
hbar_SI = 6.582119514e-16 # reduced Planck constant (eV/second)
ge_SI = 1.760859708e11 # gyromagnetic ratio of NV electron (Hz/tesla)
gc_SI = 67.28284e6 # gyromagnetic ratio of C-13 (Hz/tesla)


# physical constants in natural units (s or Hz)
alpha = 1/137.035999074 # fine structure constant
qe = np.sqrt(4*np.pi*alpha) # unit electric charge
me = 510998.928/hbar_SI # mass of electron (Hz)
ge = qe/(2*me) * 2.0023193043617 # gyromagnetic ratio of electron
gc = gc_SI * ge/ge_SI # gyromagnetic ratio of C-13

sec = 1 # one second in natural units (s)
meter = sec / c_SI # one meter in natural units (s)
nm = 1e-9*meter # one nanometer in natural units (s)

a0 = 0.35668 * nm # diamond lattice parameter (unit cell side length) at 300 K
zhat = np.array([1,1,1])/np.sqrt(3) # static magnetic field direction; NV axis

# diamond lattice vectors
ao = np.array([1,1,1])/4
a1 = np.array([0,1,1])/2
a2 = np.array([1,0,1])/2
a3 = np.array([1,1,0])/2

norm = np.linalg.norm
def hat(v): return v/norm(v)
def A(b,l,m,n):
    r = b*a0 + l*a1 + m*a2 + n*a3
    return norm( ge*gc/(4*np.pi*(norm(r)*a0)**3) * (zhat - 3*np.dot(hat(r),zhat)*hat(r)) )

# maximum allowable value of l^2, m^2, or n^2
M = (2*ge*gc / (np.pi * a0**3 * hyperfine_cutoff))**(2/3)

# list of all equivalence classes of integers (b,l,m,n)
#   i.e. sets {(b,l,m,n)} within which both l^2+m^2+n^2 and |1.5*b+l+m+n| are constant
equivalence_classes = []

# loop over each integer N < M
for N in range(1, int(M)+1):

    # sets of integers {l,m,n} satisfying l^2+m^2+n^2 == N and l^2+m^2+n^2+(l+m+n)^2 < M
    square_solution_sets = set([ (sign_l*l,sign_m*m,sign_n*n)
                                 for l in range(0, int(np.sqrt(N))+1)
                                 for m in range(l, int(np.sqrt(N))+1)
                                 for n in range(m, int(np.sqrt(N))+1)
                                 if l*l+m*m+n*n == N
                                 for sign_l in [1,-1]
                                 for sign_m in [1,-1]
                                 for sign_n in [1,-1] ])

    # if we have no solutions, skip ahead to the next value of N
    if(len(square_solution_sets) == 0): continue

    # get all unique purmutations of integer sets in square_solution_sets
    square_solutions = set(sum([ list(itertools.permutations(sss))
                                 for sss in square_solution_sets ]))

    # determine all equivalence classes of integers (b,l,m,n) for this value of N
    for b in [0,1]:
        blmn_sums = list(set([ abs(1.5*b+l+m+n) for (l,m,n) in square_solutions ]))
        for s in blmn_sums:
            equivalence_class = [ (b,l,m,n) for (l,m,n) in square_solutions
                                  if abs(1.5*b+l+m+n) == s and A(b,l,m,n) > hyperfine_cutoff ]
            if len(equivalence_class) > 1:
                equivalence_classes.append(equivalence_class)

# list of all sizes of the equivalence classes
ring_sizes = [ len(equivalence_class) for equivalence_class in equivalence_classes ]

probability = 1 - product([ (1-c13_abundance)**R * (1 + R*c13_abundance/(1-c13_abundance))
                            for R in ring_sizes ])

print(probability)
