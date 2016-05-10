#!/usr/bin/env python
import sys, os, shutil, subprocess, random, threading, time, glob

if len(sys.argv) < 7:
    print("usage: " + sys.argv[0] + " sim_type static_Bz" + \
          " c13_percentage max_cluster_size log10_samples task_num [seed]")
    exit(1)

# process inputs
sim_type = sys.argv[1]
static_Bz = sys.argv[2]
c13_percentage = sys.argv[3]
max_cluster_size = sys.argv[4]
log10_samples = sys.argv[5]
task_num = int(sys.argv[6])
assert task_num > 1
seed_text = " ".join(sys.argv[7:])

print_period = 1800 # seconds

# identify some directories and names
project_dir = os.path.dirname(os.path.realpath(__file__))
sim_file = "simulate.exe"
out_file = "data/fidelities-{}.txt".format("-".join([sim_type, static_Bz, c13_percentage,
                                                     max_cluster_size, log10_samples]))

if "TMPDIR" not in os.environ:
    print("please set the TMPDIR environment variable")
    exit(1)
job_dir = os.environ["TMPDIR"] + "/" + os.environ["USER"] + "/" + "-".join(sys.argv[1:])

# move into job directory
os.makedirs(job_dir)
shutil.copy2(project_dir+"/"+sim_file, job_dir+"/"+sim_file)
os.chdir(job_dir)

# method to execute a single simulation
def run_sample(s):
    random.seed("".join(sys.argv[1:-1]) + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(commands + seed,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    with lock:
        with open(out_file,"a") as f:
            f.write(" ".join(commands + seed)+"\n\n")
            f.write(out.decode("utf-8")+"\n")
            f.write(err.decode("utf-8")+"\n")
            f.write("-"*90 + "\n\n")

# method to copy results back into the project directory
def copy_results():
    for f in glob.glob("data/*"):
        shutil.copy2(f, project_dir+"/"+f)

# define some variavles
commands = ["./"+sim_file, "--no_output", "--"+sim_type,
            "--static_Bz", static_Bz,
            "--c13_percentage", c13_percentage,
            "--max_cluster_size", max_cluster_size]

unsigned_long_long_max = 2**64-1
print_time = time.time()
samples = int(10**float(log10_samples))
lock = threading.RLock()

# execute all simulations
for s in range(samples):

    # run each simulation in a new thread
    t = threading.Thread(target=run_sample,args=[s])
    while threading.active_count() >= task_num:
        time.sleep(1)
    t.start()

    # periodically copy results back into project directory
    if time.time() - print_time > print_period:
        copy_results()
        print_time = time.time()

# wait for all threads to finish
while threading.active_count() > 1:
    time.sleep(1)

# copy final results into the project directory
copy_results()
