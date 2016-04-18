#!/usr/bin/python3

import sys, os
import subprocess as sp
import numpy as np
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

fname = "./data/pair-test-{}-{}-{}.txt".format(start,end,log10_samples)

if os.path.exists(fname):
    cutoffs, predicted, actual = np.loadtxt(fname, unpack=True)

else:
    cutoffs = range(start,end+1)

    predicted = np.zeros(len(cutoffs))
    actual = np.zeros(len(cutoffs))
    for i in range(len(cutoffs)):
        print("starting cutoff: {} kHz".format(cutoffs[i]))
        predicted[i] = sp.check_output(["./pair-compute.py",str(cutoffs[i])])
        actual[i] = sp.check_output(["./pair-search.py",str(cutoffs[i]),
                                     str(10**log10_samples)])

    with open(fname,'w') as f:
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
plt.savefig(fname.replace(".txt",".pdf"))
