#!/usr/bin/env python
import sys, os, random, string, subprocess

if len(sys.argv) != 6:
    print("usage: " + sys.argv[0] + " sim_type static_Bz" + \
          " c13_abundance max_cluster_size log10_samples")
    exit(1)

script = os.path.abspath("fidelity-sweep.py")
args = sys.argv[1:]

job_dir = "jobs"
basename = "-".join(args)
out_file = job_dir+"/"+basename+".o"
err_file = job_dir+"/"+basename+".e"

temp_file = "".join(random.sample(string.lowercase,10))
with open(temp_file,"w") as f:
    f.write("python {} {}\n".format(script," ".join(args)))
subprocess.call(["msub", "-m", "e", "-N", basename,
                 "-o", out_file, "-e", err_file,
                 "-l", "walltime=168:00:00,pmem=100mb",
                 temp_file])
os.remove(temp_file)
