#!/usr/bin/env python
import sys, os, random, string, subprocess

if len(sys.argv) != 7:
    print("usage: " + sys.argv[0] + " sim_type static_Bz" + \
          " c13_abundance max_cluster_size log10_samples thread_cap")
    exit(1)

script = os.path.abspath("fidelity-sweep.py")
sim_args = sys.argv[1:-1]
thread_cap = sys.argv[-1]
assert int(thread_cap) > 1

job_dir = "jobs"
basename = "-".join(sim_args)
out_file = job_dir+"/"+basename+".o"
err_file = job_dir+"/"+basename+".e"

temp_file = "".join(random.sample(string.lowercase,10))
with open(temp_file,"w") as f:
    f.write("python {} {} {}\n".format(script," ".join(sim_args),thread_cap))
subprocess.call(["msub", "-m", "e", "-N", basename,
                 "-o", out_file, "-e", err_file,
                 "-l", "walltime=168:00:00,pmem=100mb,procs="+thread_cap,
                 temp_file])
os.remove(temp_file)
