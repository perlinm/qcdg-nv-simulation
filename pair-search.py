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

random.seed(str(samples)+hyperfine_cutoff)
unsigned_long_long_max = 2**64-1

with open(os.devnull,'w') as null:
    subprocess.call(["fac"],stdout=null)

    found = 0
    for s in range(samples):
        commands = ["./simulate",
                    "--pair",
                    "--hyperfine_cutoff",hyperfine_cutoff,
                    "--seed",str(random.randint(0,unsigned_long_long_max))]
        if len(sys.argv) == 4:
            commands += ["--c13_abundance",sys.argv[3]]
        pairs = subprocess.call(commands,stdout=null)
        found += 1 if pairs > 0 else 0

print(found/samples)
