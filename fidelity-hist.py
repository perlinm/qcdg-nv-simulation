#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("usage: {} fidelity_sweep_file".format(sys.argv[0]))
    exit(1)

fname = sys.argv[1]
if not os.path.isfile(fname):
    print("invalid file: {}".format(fname))
    exit(1)

results = []
collect = False
with open(fname,'r') as f:
    for line in f:
        if "--seed" in line:
            seed = int(line.split()[-1])
        if " fidelity time\n" in line:
            collect = True
            continue
        if line == "\n":
            collect = False
            continue
        if collect:
            line = line.split()
            if len(line) < 3:
                print("invalid data line")
                exit(1)
            elif len(line) == 3:
                targets = int(line[0])
            else:
                targets = [ int(target) for target in line[0:-2] ]
            results += [(seed,targets,float(line[-2]),float(line[-1]))]

results.sort(key = lambda x: x[2])
for result in results:
    print(result)

plt.hist([result[2] for result in results], color='gray')
plt.xlim(0,1)
plt.show()
