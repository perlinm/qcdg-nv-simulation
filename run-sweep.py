#!/usr/bin/env python
import sys, os, subprocess

if len(sys.argv) < 8:
    print("usage: " + sys.argv[0] + " sim_type static_Bz c13_percentage" + \
          " max_cluster_size scale_factor log10_samples walltime_in_hours [mkfac.py args]")
    exit(1)

sim_type = sys.argv[1]
static_Bz = sys.argv[2]
c13_percentage = sys.argv[3]
max_cluster_size = sys.argv[4]
scale_factor = sys.argv[5]
log10_samples = sys.argv[6]
walltime_in_hours = sys.argv[7]
mkfac_args = sys.argv[8:]
sim_args = sys.argv[1:8]

nodes = 1
tasks_per_node = 16
task_num = str(int(nodes)*tasks_per_node)

mem_per_task_estimate_in_mb = int(0.5 * 4**float(log10_samples))

project_dir = os.path.dirname(os.path.realpath(__file__))
job_dir = "jobs"
mkfac = "mkfac.py"
script = "fidelity-sweep.py"
basename = (sim_type + "-sBz-"+static_Bz + "-cp-"+c13_percentage + "-mcs-"+max_cluster_size
            + "-sf-"+scale_factor + "-ls-"+log10_samples)
job_file = job_dir+"/"+basename+".sh"
out_file = job_dir+"/"+basename+".o"

print(job_file)

resources = ["nodes={}:ppn={}".format(nodes,tasks_per_node),
             "walltime={}:00:00".format(walltime_in_hours),
             "pmem={}mb".format(mem_per_task_estimate_in_mb)]
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
