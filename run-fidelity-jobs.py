#!/usr/bin/env python
import sys, os, numpy, subprocess, glob
from basename import basename

def input_error():
    print("usage: " + sys.argv[0] + " sim_type [test] [mkfac.py args]")
    exit(1)

if len(sys.argv) < 3: input_error()

sim_type = sys.argv[1]
if sys.argv[2] == "test":
    test_jobs = True
    mkfac_args = sys.argv[3:]
else:
    test_jobs = False
    mkfac_args = sys.argv[2:]

c13_natural_percentage = 1.07

static_Bzs = [ 50, 100, 200, 500 ]
c13_factors = [ 0.1, 0.01 ]
max_cluster_sizes = [ 5 ]
scale_factors = [ 5, 8, 10 ]

project_dir = os.path.dirname(os.path.realpath(__file__))
run_script = "fidelity-sweep-job.py"

def cmd_args(sim_args, walltime):
    return ([ "{}/{}".format(project_dir,run_script) ]
            + [ str(a) for a in sim_args]
            + [ str(walltime) ] + mkfac_args)

for static_Bz in static_Bzs:
    for c13_factor in c13_factors:
        for max_cluster_size in max_cluster_sizes:
            for scale_factor in scale_factors:
                log10_samples = int(numpy.round(3 - numpy.log10(c13_factor)))
                c13_percentage = str(numpy.around(c13_natural_percentage*c13_factor,
                                                  int(3 - numpy.log10(c13_factor))))
                sim_args = [ sim_type, static_Bz, c13_percentage,
                             max_cluster_size, scale_factor, log10_samples ]

                if test_jobs:
                    subprocess.call(cmd_args(sim_args,3))

                else:
                    test_job_name = "./jobs/" + basename(sim_args) + ".o_feedback"
                    fname_candidates = glob.glob(test_job_name)
                    assert (len(fname_candidates) == 1)
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
