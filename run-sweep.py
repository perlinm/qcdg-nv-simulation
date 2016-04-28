#!/usr/bin/env python
import sys, os, random, string, subprocess

if len(sys.argv) != 8:
    print("usage: " + sys.argv[0] + " sim_type static_Bz c13_abundance" + \
          " max_cluster_size log10_samples nodes walltime")
    exit(1)

script = os.path.abspath("fidelity-sweep.py")
sim_args = sys.argv[1:6]
nodes = sys.argv[-2]
walltime = sys.argv[-1]
assert int(task_num) > 1

job_dir = "jobs"
basename = "-".join(sim_args)
out_file = job_dir+"/"+basename+".o"
err_file = job_dir+"/"+basename+".e"

tasks_per_node = 16
task_num = str(int(nodes)*tasks_per_node)
resources = "nodes={}:ppn={},walltime={},pmem=10mb".format(nodes,tasks_per_node,walltime)

temp_file = "".join(random.sample(string.lowercase,10))
with open(temp_file,"w") as f:
    f.write("python {} {} {}\n".format(script," ".join(sim_args),task_num))

subprocess.call(["msub", "-m", "e", "-N", basename,
                 "-o", out_file, "-e", err_file,
                 "-l", resources, temp_file])
os.remove(temp_file)
