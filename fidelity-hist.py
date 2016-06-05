#!/usr/bin/env python
import sys, os, re
from basename import basename

if len(sys.argv) != 2:
    print("usage: {} fidelity_data_file".format(sys.argv[0]))
    exit(1)

fname = sys.argv[1]
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

results.sort(key = lambda x: -x[0])

with open(out_file,"w") as f:
    f.write("# samples: {}\n".format(samples))
    f.write("# fidelity time seed target(s)\n")
    for result in results:
        f.write("{} {} {}".format(result[0],result[1],result[2]))
        for target in result[3]:
            f.write(" {}".format(target))
        f.write("\n")
