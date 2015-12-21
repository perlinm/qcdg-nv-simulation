#!/usr/bin/env python3

import os, sys, subprocess, random

if len(sys.argv) != 3:
    print("usage: {} samples hyperfine_cutoff (in kHz)".format(sys.argv[0]))
    exit(1)

samples = sys.argv[1]
hyperfine_cutoff = sys.argv[2]

random.seed(samples+hyperfine_cutoff)
def rnd():
    return int(random.random()*1e6)

with open(os.devnull, 'w') as null:
    subprocess.call(['fac'],stdout=null)

    found = 0
    for s in range(int(samples)):
        found += subprocess.call(["./simulate",
                                  "--pair_search",
                                  "--hyperfine_cutoff",hyperfine_cutoff,
                                  "--seed",str(rnd())],stdout=null)

print(sys.argv[2],found,samples)
