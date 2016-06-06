#!/usr/bin/env python
import sys, os, re
from basename import basename

if len(sys.argv) not in [2,3]:
    print("usage: {} fidelity_data_file [fidelity_cutoff]".format(sys.argv[0]))
    exit(1)

fname = sys.argv[1]
if len(sys.argv) == 3:
    fidelity_cutoff = float(sys.argv[2])
else:
    fidelity_cutoff = 0.5

if not os.path.isfile(fname):
    print("invalid file: {}".format(fname))
    exit(1)

out_file = "{}/fidelities-{}".format(os.path.dirname(fname),os.path.basename(fname))

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
            if len(line) == 3:
                targets = int(line[0])
            elif len(line) == 4:
                targets = [int(line[0]),int(line[1])]
            else:
                print("invalid data")
                exit(1)
            results += [(float(line[-2]),float(line[-1]),seed,targets)]

if len(results) == 0:
    print("no simulation data in file: {}".format(fname))
    exit(1)

cutoff_results = [ result for result in results if result[0] > fidelity_cutoff ]
cutoff_results.sort(key = lambda x: -x[0])

with open(out_file,"w") as f:
    f.write("# samples: {}\n".format(samples))
    f.write("# total results: {}\n".format(len(results)))
    f.write("# fidelity time seed target(s)\n")
    for result in cutoff_results:
        f.write("{} {} {}".format(result[0],result[1],result[2]))
        for target in result[3]:
            f.write(" {}".format(target))
        f.write("\n")
