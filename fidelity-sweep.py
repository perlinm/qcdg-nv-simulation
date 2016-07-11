#!/usr/bin/env python
import sys, os, shutil, subprocess, random, threading, time, glob
from basename import basename

if len(sys.argv) < 8:
    print("usage: " + sys.argv[0] + " [sim_args...] task_num")
    exit(1)

# process inputs
task_num = int(sys.argv[-1])
assert task_num > 1
del sys.argv[-1]

log10_samples_hook = "--log10_samples"
log10_samples = sys.argv[sys.argv.index(log10_samples_hook)+1]
del sys.argv[sys.argv.index(log10_samples_hook)+1]
del sys.argv[sys.argv.index(log10_samples_hook)]

sim_args = sys.argv[1:]

print_period = 1800 # seconds

# identify some directories and names
project_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = "data"
sim_file = "simulate.exe"
summary_script = "fidelity-summary.py"
out_file = "{}/{}.txt".format(data_dir,basename(sim_args))
base_cmd = ["./"+sim_file] + sim_args

if "TMPDIR" not in os.environ:
    print("please set the TMPDIR environment variable")
    exit(1)
job_dir = os.environ["TMPDIR"] + "/" + os.environ["USER"] + "/" + "-".join(sys.argv[1:])

# move into job directory
os.makedirs(job_dir)
shutil.copy2(project_dir+"/"+sim_file, job_dir+"/"+sim_file)
os.chdir(job_dir)
os.mkdir(data_dir)

# method to execute a single simulation
def run_sample(s):
    random.seed("".join(sim_args) + str(s))
    seed = ["--seed", str(random.randint(0,unsigned_long_long_max))]
    process = subprocess.Popen(base_cmd + seed,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    with lock:
        with open(out_file,"a") as f:
            f.write(" ".join(base_cmd + seed)+"\n\n")
            f.write(out.decode("utf-8")+"\n")
            f.write(err.decode("utf-8")+"\n")
            f.write("-"*90 + "\n\n")

# method to copy results back into the project directory
def copy_results():
    for f in glob.glob(data_dir+"/*"):
        shutil.copy2(f, project_dir+"/"+f)

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

# generate fidelity summary file
os.chdir(project_dir)
subprocess.call(["./"+summary_script,out_file])

