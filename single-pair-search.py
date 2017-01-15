#!/usr/bin/env python
import sys, os, subprocess, random, threading

if len(sys.argv) != 8:
    print("usage: " + sys.argv[0] + " hyperfine_cutoff_in_kHz c13_factor" + \
          " min_hyperfine_xy_in_kHz nuclear_isolation_in_Hz" + \
          " larmor_isolation_in_kHz samples task_num")
    exit(1)

hyperfine_cutoff = sys.argv[1]
c13_factor = sys.argv[2]
min_hyperfine_xy = sys.argv[3]
nuclear_isolation_in_Hz = sys.argv[4]
larmor_isolation = sys.argv[5]
samples = int(sys.argv[6])
task_num = int(sys.argv[7])
assert task_num >= 1

sim_file = "{}/simulate.exe".format(os.path.dirname(os.path.realpath(__file__)))

unsigned_long_long_max = 2**64-1

commands = [sim_file, "--pair_search",
            "--c13_factor", c13_factor,
            "--hyperfine_cutoff", hyperfine_cutoff,
            "--min_hyperfine_xy", min_hyperfine_xy,
            "--nuclear_isolation", nuclear_isolation_in_Hz,
            "--larmor_isolation", larmor_isolation]

lock = threading.RLock()

samples_with_pairs = 0
hyperfine_xy_values = []

def run_sample(s):
    random.seed("".join(sys.argv[1:-1]) + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed, stdout=subprocess.PIPE)
    out,_ = process.communicate()
    out = out.decode("utf-8")
    pairs = []
    start_reading_pairs = False
    for line in out.split("\n"):
        if "idx1 idx2 hyperfine_xy" in line:
            start_reading_pairs = True
            continue
        if start_reading_pairs and len(line.split()) == 3:
            pairs.append(line.split())

    with lock:
        global samples_with_pairs, hyperfine_xy_values
        hyperfines = [ float(pair[2]) for pair in pairs ]
        if len(hyperfines) > 0:
            samples_with_pairs += 1

for s in range(samples):
    t = threading.Thread(target=run_sample,args=[s])
    while threading.active_count() > task_num: None
    t.start()

while threading.active_count() > 1: None

print(samples_with_pairs/samples)
