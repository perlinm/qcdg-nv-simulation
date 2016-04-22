#!/usr/bin/env python3
import sys, os, subprocess, random

if len(sys.argv) != 4:
    print("usage: {} hyperfine_cutoff_in_kHz samples c13_abundance".format(sys.argv[0]))
    exit(1)

hyperfine_cutoff = sys.argv[1]
samples = int(sys.argv[2])
c13_abundance = float(sys.argv[3])

work_dir = os.path.dirname(os.path.realpath(__file__))
sim_name = "simulate"
sim_file = work_dir + "/" + sim_name

subprocess.call(["fac"])

random.seed(hyperfine_cutoff + str(samples) + str(c13_abundance))
unsigned_long_long_max = 2**64-1

found = 0
commands = [sim_file, "--no_output", "--pair",
            "--hyperfine_cutoff", str(hyperfine_cutoff),
            "--c13_abundance", str(c13_abundance)]
for s in range(samples):
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed, stdout=subprocess.PIPE)
    process.communicate()
    pairs = process.returncode
    found += 1 if pairs > 0 else 0

print(found/samples)
