#!/usr/bin/env python
import sys, os, subprocess, numpy
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    print("usage: {} cutoff_start cutoff_end log10_samples".format(sys.argv[0]))
    exit(1)

start = int(sys.argv[1])
end = int(sys.argv[2])
log10_samples = int(sys.argv[3])

if not start < end:
    print("cutoff start must be less than end")
    exit(2)

work_dir = os.path.dirname(os.path.realpath(__file__))
out_name = "data/pairs-{}-{}-{}.txt".format(start,end,log10_samples)
out_file = work_dir + "/" + out_name

if os.path.exists(out_file):
    cutoffs, predicted, actual = numpy.loadtxt(out_file, unpack=True)

else:
    cutoffs = range(start,end+1)

    predicted = numpy.zeros(len(cutoffs))
    actual = numpy.zeros(len(cutoffs))
    for i in range(len(cutoffs)):
        print("starting cutoff: {} kHz".format(cutoffs[i]))
        predicted[i] = subprocess.check_output(["./pair-compute.py",str(cutoffs[i])])
        actual[i] = subprocess.check_output(["./pair-search.py",str(cutoffs[i]),
                                             str(10**log10_samples)])

    with open(out_file,'w') as f:
        f.write("# log10_samples: {}\n".format(log10_samples))
        f.write("# hyperfine_cutoff predicted actual\n")
        for i in range(len(cutoffs)):
            f.write("{} {} {}\n".format(cutoffs[i],predicted[i],actual[i]))

plt.title("Larmor pair probability test with $10^{{{}}}$ samples".format(log10_samples))
plt.plot(cutoffs,predicted,"k-",label="predicted")
plt.plot(cutoffs,actual,"k.",label="found")
plt.xlabel("Hyperfine cutoff [kHz]")
plt.ylabel("Proportion")
plt.ylim(0,1)
plt.legend(loc="best")
plt.savefig(out_file.replace(".txt",".pdf"))
