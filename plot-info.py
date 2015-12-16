#!/usr/bin/python3
from pylab import *

index, w, A, A_perp, dw_min, f_DD, operation_time, fidelity = \
       loadtxt('iswap_info.txt', unpack=True)

mxi = fidelity == max(fidelity)

print(fidelity[mxi])
print(operation_time[mxi])

figure()
title('fidelity histogram (scale$=100$)')
hist(fidelity[fidelity > 0.9])
xlabel('fidelity')

show()
