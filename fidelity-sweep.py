#!/usr/bin/env python3
import sys, os, random
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 6:
    print("usage: " + sys.argv[0] + " sim_type static_Bz" + \
          " c13_abundance max_cluster_size log10_samples [seed]")
    exit(1)

sim_type = sys.argv[1]
static_Bz = int(sys.argv[2])
c13_abundance = float(sys.argv[3])
max_cluster_size = int(sys.argv[4])
log10_samples = int(sys.argv[5])
seed_text = ' '.join(sys.argv[5:])

fname = "./data/fidelities-{}-{}-{}-{}-{}.txt".format(sim_type, static_Bz, c13_abundance,
                                                      max_cluster_size, log10_samples)

random.seed(fname + seed_text)
unsigned_long_long_max = 2**64-1

samples = int(10**log10_samples)
sp.call(["fac"])
with open(fname,'w') as output:
    commands = ["./simulate", "--no_output", "--" + sim_type,
                "--static_Bz", str(static_Bz),
                "--c13_abundance", str(c13_abundance),
                "--max_cluster_size", str(max_cluster_size)]
    for s in range(samples):
        print("{} / {}".format(s,samples))
        seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
        output.write(' '.join(commands + seed) + "\n\n")
        output.flush()
        process = sp.Popen(commands + seed, stdout=sp.PIPE)
        out, err = process.communicate()
        output.write(out.decode("utf-8") + "\n")


