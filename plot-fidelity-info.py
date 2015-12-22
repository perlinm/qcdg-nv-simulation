#!/usr/bin/python3
from pylab import *

if len(sys.argv) != 2:
    print("useage: %s file" % sys.argv[0])
    exit(1)

fname = sys.argv[1]

if "fidelities" not in fname:
    print("invalid input file")
    exit(1)

index, w, A, A_perp, dw_min, f_DD, operation_time, fidelity = loadtxt(fname, unpack=True)

mxi = fidelity == max(fidelity)
print("max fidelity:",max(fidelity))
print("operation time(s):",operation_time[mxi])

cutoff = 0.95
print("fidelities above {}: {}".format(cutoff,fidelity[fidelity > cutoff]))

f2 = multiply(fidelity,fidelity)

figure()
hist(f2[f2 > 0.5])

show()
