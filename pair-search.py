#!/usr/bin/env python
import sys, os, subprocess, random, threading

if len(sys.argv) != 5:
    print("usage: " + sys.argv[0] + " hyperfine_cutoff_in_kHz" + \
          " samples c13_percentage thread_cap")
    exit(1)

hyperfine_cutoff = sys.argv[1]
samples = int(sys.argv[2])
c13_percentage = float(sys.argv[3])
thread_cap = int(sys.argv[4])
assert thread_cap > 1

work_dir = os.path.dirname(os.path.realpath(__file__))
sim_name = "simulate.exe"
sim_file = work_dir + "/" + sim_name

unsigned_long_long_max = 2**64-1

commands = [sim_file, "--no_output", "--pair",
            "--hyperfine_cutoff", str(hyperfine_cutoff),
            "--c13_percentage", str(c13_percentage)]

lock = threading.RLock()

found = 0
def run_sample(s):
    random.seed("".join(sys.argv[1:-1]) + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed, stdout=subprocess.PIPE)
    process.communicate()
    pairs = process.returncode
    with lock:
        global found
        found += 1 if pairs > 0 else 0

for s in range(samples):
    t = threading.Thread(target=run_sample,args=[s])
    while threading.active_count() >= thread_cap: None
    t.start()

print(found/samples)
