#!/usr/bin/env python
import sys, os, subprocess, random, threading, time

if len(sys.argv) < 7:
    print("usage: " + sys.argv[0] + " sim_type static_Bz" + \
          " c13_abundance max_cluster_size log10_samples thread_cap [seed]")
    exit(1)

sim_type = sys.argv[1]
static_Bz = int(sys.argv[2])
c13_abundance = float(sys.argv[3])
max_cluster_size = int(sys.argv[4])
log10_samples = int(sys.argv[5])
thread_cap = int(sys.argv[6])
assert thread_cap > 1
seed_text = ' '.join(sys.argv[7:])

work_dir = os.path.dirname(os.path.realpath(__file__))
out_name = "data/fidelities-{}-{}-{}-{}-{}.txt".format(sim_type, static_Bz, c13_abundance,
                                                       max_cluster_size, log10_samples)
sim_name = "simulate.exe"
out_file = work_dir + "/" + out_name
sim_file = work_dir + "/" + sim_name

unsigned_long_long_max = 2**64-1
samples = int(10**log10_samples)

commands = [sim_file, "--no_output", "--" + sim_type,
            "--static_Bz", str(static_Bz),
            "--c13_abundance", str(c13_abundance),
            "--max_cluster_size", str(max_cluster_size)]

lock = threading.RLock()

def run_sample(s):
    random.seed("".join(sys.argv[1:-1]) + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    with lock:
        with open (out_file,"a") as f:
            f.write(' '.join(commands + seed)+"\n\n")
            f.write(out.decode("utf-8")+"\n")
            f.write(err.decode("utf-8")+"\n")
            f.write("-"*80 + "\n\n")

for s in range(samples):
    print("{} / {}".format(s,samples))
    t = threading.Thread(target=run_sample,args=[s])
    while threading.active_count() >= thread_cap:
        time.sleep(1)
    t.start()
