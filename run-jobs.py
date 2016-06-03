#!/usr/bin/env python
import sys, numpy, subprocess, glob
from basename import basename

test_jobs = "test" in sys.argv
real_jobs = "real" in sys.argv
assert (test_jobs or real_jobs)

sim_type = "swap_nvst"
c13_base_percentage = 1.07

def cmd_args(sim_args, walltime):
    return [ "./run-sweep.py" ] + [ str(a) for a in sim_args] + [ walltime, "whide" ]

for static_Bz in [ 500, 1000, 1500 ]:
    for c13_factor in [ 1, 0.1, 0.01 ]:
        for max_cluster_size in [ 4 ]:
            for scale_factor in [ 5, 10, 20 ]:
                log10_samples = int(numpy.round(3 - numpy.log10(c13_factor)))
                c13_percentage = str(numpy.around(c13_base_percentage*c13_factor,
                                                  int(3 - numpy.log10(c13_factor))))
                sim_args = [ sim_type, static_Bz, c13_percentage,
                             max_cluster_size, scale_factor, log10_samples ]

                if test_jobs:
                    subprocess.call(cmd_args(sim_args,3))

                if real_jobs:
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
