#!/usr/bin/env python
import sys, os, tempfile, subprocess

if len(sys.argv) != 7:
    print("usage: " + sys.argv[0] + " sim_type static_Bz c13_abundance" + \
          " max_cluster_size log10_samples walltime_in_days")
    exit(1)

script = os.path.abspath("fidelity-sweep.py")
sim_args = sys.argv[1:6]
walltime_in_days = sys.argv[-1]

job_dir = "jobs"
basename = "-".join(sim_args)
out_file = job_dir+"/"+basename+".o"
err_file = job_dir+"/"+basename+".e"

nodes = 1
tasks_per_node = 16
task_num = str(int(nodes)*tasks_per_node)
resources = ",".join(["nodes={}:ppn={}".format(nodes,tasks_per_node),
                      "walltime={}:00:00:00".format(walltime_in_days),
                      "pmem=1mb"])

_,temp_file = tempfile.mkstemp()
with open(temp_file,"w") as f:
    f.write("python {} {} {}\n".format(script," ".join(sim_args),task_num))

os.system("cat %s"%temp_file)
subprocess.call(["msub", "-m", "e", "-N", basename,
                 "-o", out_file, "-e", err_file,
                 "-l", resources, temp_file])
os.remove(temp_file)
