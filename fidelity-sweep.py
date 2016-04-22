#!/usr/bin/env python3
import sys, os, subprocess, random

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

work_dir = os.path.dirname(os.path.realpath(__file__))
out_name = "data/fidelities-{}-{}-{}-{}-{}.txt".format(sim_type, static_Bz, c13_abundance,
                                                       max_cluster_size, log10_samples)
sim_name = "simulate"

out_file = work_dir + "/" + out_name
sim_file = work_dir + "/" + sim_name

subprocess.call(["fac"])

random.seed(out_name + seed_text)
unsigned_long_long_max = 2**64-1

samples = int(10**log10_samples)
with open(out_file,'w') as output:
    commands = [sim_file, "--no_output", "--" + sim_type,
                "--static_Bz", str(static_Bz),
                "--c13_abundance", str(c13_abundance),
                "--max_cluster_size", str(max_cluster_size)]
    for s in range(samples):
        print("{} / {}".format(s,samples))
        seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
        output.write(' '.join(commands + seed) + "\n\n")
        output.flush()
        process = subprocess.Popen(commands + seed,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        if process.returncode == 0:
            output.write(out.decode("utf-8") + "\n")


