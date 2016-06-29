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
        if not collect:
            if "--seed" in line:
                seed = int(line.split()[line.split().index("--seed")+1])
            elif " fidelity" in line:
                description = line.replace("\n"," seed")
                collect = True
        else:
            if line == "\n":
                collect = False
                continue
            results += [ "{} {}".format(line.replace("\n",""),seed) ]

if len(results) == 0:
    print("no simulation data in file: {}".format(fname))
    exit(1)

fidelity_column = description.split().index("fidelity")
cutoff_results = [ result for result in results
                   if float(result.split()[fidelity_column]) > fidelity_cutoff ]
cutoff_results.sort(key = lambda r: -float(r.split()[fidelity_column]))

with open(out_file,"w") as f:
    f.write("# samples: {}\n".format(samples))
    f.write("# total results: {}\n".format(len(results)))
    f.write("# {}\n".format(description))
    for result in cutoff_results:
        f.write(result+"\n")
