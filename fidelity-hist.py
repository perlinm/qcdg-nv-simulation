#!/usr/bin/env python
import sys, os, re

if len(sys.argv) not in [2,3]:
    print("usage: {} fidelity_data_file [cutoff] [print]".format(sys.argv[0]))
    exit(1)

print_results = True if sys.argv[-1] == "print" else False
if print_results: del sys.argv[-1]

fname = sys.argv[1]
if len(sys.argv) == 3:
    cutoff = float(sys.argv[2])
else:
    cutoff = 0.9

if not os.path.isfile(fname):
    print("invalid file: {}".format(fname))
    exit(1)
samples = 10**int(re.split("-|\.",fname)[-2])

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
            results += [(float(line[-2]),float(line[-1]),targets,seed)]

results.sort(key = lambda x: x[0])
filtered_results = [ r for r in results if r[0] > cutoff ]

if print_results:
    for result in filtered_results:
        print(result)

print(len(filtered_results)/len(results))
print(len(filtered_results)/samples)
