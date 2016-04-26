#!/usr/bin/env python
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
seed_text = ' '.join(sys.argv[6:])

work_dir = os.path.dirname(os.path.realpath(__file__))
out_name = "data/fidelities-{}-{}-{}-{}-{}.txt".format(sim_type, static_Bz, c13_abundance,
                                                       max_cluster_size, log10_samples)
sim_name = "simulate.exe"
out_file = work_dir + "/" + out_name
sim_file = work_dir + "/" + sim_name

unsigned_long_long_max = 2**64-1
samples = int(10**log10_samples)

commands = [sim_file, "--no_output", "--" + sim_type,
            "--static_Bz", str(static_Bz),
            "--c13_abundance", str(c13_abundance),
            "--max_cluster_size", str(max_cluster_size)]

def run_sample(s):
    random.seed(out_name + seed_text + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    output_text = ' '.join(commands + seed) + "\n\n"
    output_text += out.decode("utf-8")+"\n"
    output_text += err.decode("utf-8")+"\n"
    output_text += "----------------------------------------------------------------------\n\n"
    return output_text

with open(out_file,'w') as output:
    for s in range(samples):
        print("{} / {}".format(s,samples))
        text = run_sample(s)
        output.write(text)
        output.flush()

