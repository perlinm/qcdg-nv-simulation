#!/usr/bin/env python3

import os, sys, subprocess, random

if len(sys.argv) not in [2,3,4]:
    print("usage: {} hyperfine_cutoff_in_kHz [samples] [c13_abundance]".format(sys.argv[0]))
    exit(1)

hyperfine_cutoff = sys.argv[1]
try:
    samples = int(sys.argv[2])
except:
    samples = 100
try:
    c13_abundance = float(sys.argv[3])
except:
    c13_abundance = 1.07

random.seed(hyperfine_cutoff + str(samples) + str(c13_abundance))
unsigned_long_long_max = 2**64-1

with open(os.devnull,'w') as null:
    subprocess.call(["fac"], stdout=null)

    commands = ["./simulate", "--no_output", "--pair",
                "--hyperfine_cutoff", str(hyperfine_cutoff),
                "--c13_abundance", str(c13_abundance)]
    found = 0
    for s in range(samples):
        seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
        pairs = subprocess.call(commands + seed, stdout=null)
        found += 1 if pairs > 0 else 0

print(found/samples)
