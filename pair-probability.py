#!/usr/bin/env python
import sys, os, subprocess, numpy

if len(sys.argv) not in [5,6,7]:
    print("usage: " + sys.argv[0] + " cutoff_start cutoff_end" + \
          " c13_percentage log10_samples [task_num] [plot]")
    exit(1)

make_plot = "plot" in sys.argv
if make_plot:
    sys.argv.remove("plot")
    import matplotlib.pyplot as plt

start = int(sys.argv[1])
end = int(sys.argv[2])
c13_percentage = float(sys.argv[3])
log10_samples = int(sys.argv[4])
try:
    task_num = int(sys.argv[5])
except:
    task_num = 2
assert task_num > 1

if not start < end:
    print("cutoff_start must be less than cutoff_end")
    exit(2)

project_dir = os.path.dirname(os.path.realpath(__file__))
out_name = "pairs-{}-{}-{}.txt".format(start,end,log10_samples)
out_file = project_dir + "/data/" + out_name

if os.path.exists(out_file):
    cutoffs, predicted, actual = numpy.loadtxt(out_file, unpack=True)

else:
    cutoffs = range(start,end+1)

    predicted = numpy.zeros(len(cutoffs))
    actual = numpy.zeros(len(cutoffs))
    for i in range(len(cutoffs)):
        print("starting cutoff: {} kHz".format(cutoffs[i]))
        compute_cmds = ["./pair-compute.py",str(cutoffs[i]),str(c13_percentage)]
        search_cmds = ["./pair-search.py",str(cutoffs[i]),str(c13_percentage),
                       str(10**log10_samples),str(task_num)]
        predicted[i] = subprocess.check_output(compute_cmds)
        actual[i] = subprocess.check_output(search_cmds)

    with open(out_file,'w') as f:
        f.write("# log10_samples: {}\n".format(log10_samples))
        f.write("# hyperfine_cutoff predicted actual\n")
        for i in range(len(cutoffs)):
            f.write("{} {} {}\n".format(cutoffs[i],predicted[i],actual[i]))

if make_plot:
    plt.title("Larmor pair probability test" + \
              " with $10^{{{}}}$ samples".format(log10_samples))
    plt.plot(cutoffs,predicted,"k-",label="predicted")
    plt.plot(cutoffs,actual,"k.",label="found")
    plt.xlabel("Hyperfine cutoff [kHz]")
    plt.ylabel("Proportion")
    plt.ylim(0,1)
    plt.legend(loc="best")
    plt.savefig(out_file.replace(".txt",".pdf"))
    plt.show()
