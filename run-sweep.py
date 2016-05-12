#!/usr/bin/env python
import sys, os, shutil, tempfile, subprocess

if len(sys.argv) < 7:
    print("usage: " + sys.argv[0] + " sim_type static_Bz c13_abundance" + \
          " max_cluster_size log10_samples walltime_in_days [mkfac.py args]")
    exit(1)

sim_args = sys.argv[1:6]
walltime_in_days = sys.argv[6]
mkfac_args = sys.argv[7:]

project_dir = os.path.dirname(os.path.realpath(__file__))
job_dir = "jobs"
mkfac = "mkfac.py"
script = "fidelity-sweep.py"
basename = "-".join(sim_args)
job_file = job_dir+"/"+basename+".sh"
out_file = job_dir+"/"+basename+".o"

nodes = 1
tasks_per_node = 16
task_num = str(int(nodes)*tasks_per_node)

resources = ["nodes={}:ppn={}".format(nodes,tasks_per_node),
             "walltime={}:00:00:00".format(walltime_in_days),
             "pmem=50mb"]
options = ["-N "+basename, "-m e", "-j oe", "-o "+out_file]
for resource in resources:
    options += ["-l "+resource]

job_text = "#!/usr/bin/env sh\n"
for option in options:
    job_text += "#MSUB {}\n".format(option)
job_text += "\n"
job_text += "python {}/{} {} {}\n".format(project_dir,script," ".join(sim_args),task_num)

with open(job_file,"w") as f: f.write(job_text)

subprocess.call(["./"+mkfac] + mkfac_args)
subprocess.call(["msub",job_file])
