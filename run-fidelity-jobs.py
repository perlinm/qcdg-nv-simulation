#!/usr/bin/env python
import sys, os, numpy, subprocess, glob
from basename import *

test_flag = "test"
whide_flag = "whide"

if len(sys.argv) < 2:
    print("usage: {} [{}] [{}] --sim_type [sim_opts...]"
          .format(sys.argv[0], test_flag,whide_flag))
    exit(1)

test_jobs = test_flag in sys.argv
if test_jobs: sys.argv.remove(test_flag)
whide = whide_flag in sys.argv
if whide: sys.argv.remove(whide_flag)

sim_type = sys.argv[1]
sim_opts = sys.argv[2:]

static_Bzs = [ 50, 100, 200, 500 ]
c13_factors = [ 0.1, 0.01 ]
max_cluster_sizes = [ 5 ]
scale_factors = [ 10, 13, 15 ]

project_dir = os.path.dirname(os.path.realpath(__file__))
run_script = "fidelity-sweep-job.py"

def cmd_args(sim_args, walltime):
    return ([ "{}/{}".format(project_dir,run_script) ]
            + ([ whide_flag ] if whide else [])
            + [ str(walltime) ]
            + [ str(a) for a in sim_args])

for static_Bz in static_Bzs:
    for c13_factor in c13_factors:
        for max_cluster_size in max_cluster_sizes:
            for scale_factor in scale_factors:
                log10_samples = int(numpy.round(3 - numpy.log10(c13_factor)))

                sim_args = [ sim_type,
                             "--static_Bz", static_Bz,
                             "--c13_factor", c13_factor,
                             "--max_cluster_size", max_cluster_size,
                             "--scale_factor", scale_factor ] + sim_opts \
                             + [ log10_samples_hook, log10_samples ]

                if test_jobs:
                    subprocess.call(cmd_args(sim_args,3))

                else:
                    test_job_name = "./jobs/" + basename(sim_args) + ".o_feedback"
                    fname_candidates = glob.glob(test_job_name)
                    if not len(fname_candidates) == 1:
                        print("please run test jobs first")
                        exit(1)
                    fname = fname_candidates[0]

                    with open(fname,"r") as f:
                        for line in f:
                            if "WallTime" in line:
                                test_time_text = line.split()[2]
                                break

                    hh_mm_ss = [ int(t) for t in test_time_text.split(":") ]
                    test_time = hh_mm_ss[-1] + hh_mm_ss[-2]*60 + hh_mm_ss[-3]*3600

                    time_factor = 2
                    sim_args[-1] += time_factor
                    job_time = int(numpy.ceil((10**time_factor * 1.5 * test_time)/3600))
                    job_dd_hh = str(job_time//24) + ":" + str(job_time%24).zfill(2)

                    subprocess.call(cmd_args(sim_args,job_dd_hh))
