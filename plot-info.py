#!/usr/bin/python3
from pylab import *

index, w, A, A_perp, dw_min, f_DD, operation_time, fidelity = \
       loadtxt('iswap_info.txt', unpack=True)

mxi = fidelity == max(fidelity)

print('max fidelity:',max(fidelity))
print('operation time(s):',operation_time[mxi])
print()
print(fidelity[fidelity > 0.95])

f2 = multiply(fidelity,fidelity)

figure()
hist(f2[f2 > 0.5])

show()
