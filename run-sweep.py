#!/usr/bin/env python
import sys, os, shutil, tempfile, subprocess

if len(sys.argv) != 7:
    print("usage: " + sys.argv[0] + " sim_type static_Bz c13_abundance" + \
          " max_cluster_size log10_samples walltime_in_days")
    exit(1)

sim_args = sys.argv[1:6]
walltime_in_days = sys.argv[-1]

work_dir = os.path.dirname(os.path.realpath(__file__))
script = work_dir + "/fidelity-sweep.py"
basename = "-".join(sim_args)
out_file = "jobs/"+basename+".o"

nodes = 1
tasks_per_node = 16
task_num = str(int(nodes)*tasks_per_node)
resources = ",".join(["nodes={}:ppn={}".format(nodes,tasks_per_node),
                      "walltime={}:00:00:00".format(walltime_in_days),
                      "pmem=1mb"])

_,temp_file = tempfile.mkstemp()
with open(temp_file,"w") as f:
    f.write("python {} {} {}\n".format(script," ".join(sim_args),task_num))

subprocess.call(["msub", "-m", "e", "-N", basename,
                 "-j", "oe", "-o", out_file,
                 "-l", resources, temp_file])
os.remove(temp_file)
