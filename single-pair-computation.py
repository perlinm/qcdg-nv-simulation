#!/usr/bin/env python
import sys
import numpy as np

if len(sys.argv) != 3:
    print("usage: {} hyperfine_cutoff_in_kHz c13_factor".format(sys.argv[0]))
    exit(1)

# process inputs
hyperfine_cutoff = float(sys.argv[1])*1000  # cutoff for hyperfine field strength; Hz
c13_abundance = float(sys.argv[2])*0.0107

# physical constants in SI
c_SI = 299792458. # speed of light (meters/second)
hbar_SI = 6.582119514e-16 # reduced Planck constant (eV/second)
ge_SI = -1.760859708e11 # gyromagnetic ratio of NV electron (Hz/tesla)
gC13_SI = 67.28284e6 # gyromagnetic ratio of C-13 (Hz/tesla)

# physical constants in natural units (s or Hz)
alpha = 1/137.035999074 # fine structure constant
qe = np.sqrt(4*np.pi*alpha) # unit electric charge
me = 510998.928/hbar_SI # mass of electron (Hz)
ge = -qe/(2*me) * 2.0023193043617 # gyromagnetic ratio of electron
gC13 = gC13_SI * ge/ge_SI # gyromagnetic ratio of C-13

sec = 1. # one second in natural units (s)
meter = sec / c_SI # one meter in natural units (s)
nm = 1e-9*meter # one nanometer in natural units (s)

a0 = 0.35668 * nm # diamond lattice parameter (unit cell side length) at 300 K
zhat = np.array([1,1,1])/np.sqrt(3) # static magnetic field direction; NV axis

# diamond lattice vectors
ao = np.array([1,1,1])/2.
a1 = np.array([0,1,1])
a2 = np.array([1,0,1])
a3 = np.array([1,1,0])

norm = np.linalg.norm
def hat(v): return v/norm(v)
def A(b,l,m,n):
    r = b*ao + l*a1 + m*a2 + n*a3
    return norm( ge*gC13/(4*np.pi*(norm(r)*a0/2)**3) *
                 (zhat-3*np.dot(hat(r),zhat)*hat(r)) )

def z_projection_sum(iv): return abs(3*iv[0] + 4*(iv[1]+iv[2]+iv[3]))
def xy_projection_sum(iv): return 3*(iv[1]**2+iv[2]**2+iv[3]**2) - (iv[1]+iv[2]+iv[3])**2

# maximum allowable magnitude of l, m, or n
max_lmn = int((2*abs(ge*gC13) / (np.pi * a0**3 * hyperfine_cutoff))**(1./3)+1./2)

# list of all larmor sets
larmor_sets = []

# loop over each integer z_sum up to the maximum possible value given the hyperfine cutoff
# start at 1 to exclude nuclei which lie in the x-y plane
for z_sum in range(1,12*max_lmn+1):

    # vector of integers (b,l,m,n) with this z_sum
    z_sum_solutions = set([ (b,l,m,n)
                            for b in [0,1]
                            for l in range(-max_lmn,max_lmn+1)
                            for m in range(-max_lmn,max_lmn+1)
                            for n in range(-max_lmn,max_lmn+1)
                            if z_projection_sum((b,l,m,n)) == z_sum
                            if A(b,l,m,n) > hyperfine_cutoff ])

    # determine all xy_sums for these solutions of z_sum
    xy_sums = set([ xy_projection_sum(z_solution) for z_solution in z_sum_solutions ])

    # collect solutions of z_sum into larmor sets if they have the same xy_sum
    for xy_sum in xy_sums:
        larmor_set = [ z_solution for z_solution in z_sum_solutions
                       if xy_projection_sum(z_solution) == xy_sum ]
        larmor_sets.append(larmor_set)

# return a the sizes of the nontrivial parallel subsets in a given equivalence class
def parallel_subset_sizes(larmor_set):
    r_xy_vecs = [ np.array([-2*l+m+n,-2*m+n+l,-2*n+l+m]) for _,l,m,n in larmor_set ]

    parallel_sets = []
    added_to_set = []
    for i in range(len(r_xy_vecs)):
        if i in added_to_set: continue
        parallel_sets += [[i]]
        added_to_set += [i]

        r_xy = r_xy_vecs[i]
        for j in range(i+1,len(r_xy_vecs)):
            if j in added_to_set: continue

            s_xy = r_xy_vecs[j]
            if np.array_equal(r_xy,s_xy) or np.array_equal(r_xy,-s_xy):
                parallel_sets[-1] += [j]
                added_to_set += [j]
    parallel_set_sizes = [ len(parallel_set) for parallel_set in parallel_sets ]
    return [ size for size in parallel_set_sizes if size > 1 ]

larmor_set_info = [ (len(larmor_set), parallel_subset_sizes(larmor_set))
                    for larmor_set in larmor_sets ]

# return product of elements in a list
def product(list):
    out = 1
    for el in list: out *= el
    return out

# choose function (mutiplicative formula)
def nCk(n,k): return product([ n+1-i for i in range(1,k+1)]) / product(range(1,k+1))

# probability of having at least one larmor set with exactly two addressable nuclei
probability = 1 - product([ 1 - c13_abundance**2 * (1-c13_abundance)**(R-2) *
                            ( nCk(R,2) - sum([ nCk(s,2) for s in pss ]) )
                            for eci in larmor_set_info
                            for R in [eci[0]]
                            for pss in [eci[1]] ])
print(probability)
