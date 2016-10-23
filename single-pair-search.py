#!/usr/bin/env python
import sys, os, subprocess, random, threading

if len(sys.argv) != 6:
    print("usage: " + sys.argv[0] + " hyperfine_cutoff_in_kHz" + \
          " c13_factor pair_isolation samples task_num")
    exit(1)

hyperfine_cutoff = sys.argv[1]
c13_factor = float(sys.argv[2])
pair_isolation = int(sys.argv[3])
samples = int(sys.argv[4])
task_num = int(sys.argv[5])
assert task_num >= 1

sim_file = "{}/simulate.exe".format(os.path.dirname(os.path.realpath(__file__)))

unsigned_long_long_max = 2**64-1

commands = [sim_file, "--pair_search",
            "--hyperfine_cutoff", str(hyperfine_cutoff),
            "--c13_factor", str(c13_factor),
            "--pair_isolation", str(pair_isolation)]

lock = threading.RLock()

samples_with_pairs = 0
hyperfine_perp_values = []

def run_sample(s):
    random.seed("".join(sys.argv[1:-1]) + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed, stdout=subprocess.PIPE)
    out,_ = process.communicate()
    out = out.decode("utf-8")
    pairs = []
    start_reading_pairs = False
    for line in out.split("\n"):
        if "idx1 idx2 hyperfine_perp" in line:
            start_reading_pairs = True
            continue
        if start_reading_pairs and len(line.split()) == 3:
            pairs.append(line.split())

    with lock:
        global samples_with_pairs, hyperfine_perp_values
        samples_with_pairs += 1 if len(pairs) > 0 else 0
        hyperfine_perp_values += [ float(pair[2]) for pair in pairs ]

for s in range(samples):
    t = threading.Thread(target=run_sample,args=[s])
    while threading.active_count() > task_num: None
    t.start()

while threading.active_count() > 1: None

if sum(hyperfine_perp_values) == 0:
    mean_hyperfine_perp = 0
else:
    mean_hyperfine_perp = sum(hyperfine_perp_values)/samples_with_pairs

print(samples_with_pairs/samples, "%.3g" % mean_hyperfine_perp)
