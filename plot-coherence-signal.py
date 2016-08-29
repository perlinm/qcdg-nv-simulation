#!/usr/bin/env python
import sys, matplotlib
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) not in [2,3]:
    print('useage: %s file [show]' % sys.argv[0])
    exit(1)

fname = sys.argv[1]

show = (len(sys.argv) == 3)

# return number of valid floats in a line of text
def nums_in(line):
    n = 0
    for l in line.split():
        try:
            float(l.strip())
            n += 1
        except: None
    return n

def make_plot(target,k_DD,f_DD,coherence):

    plt.plot(f_DD,coherence,'k-',linewidth=2)

    plt.xlim(0,f_DD[-1])
    plt.ylim(-1,1)

    plt.xlabel("$f_{}$".format(k_DD))
    plt.ylabel("Coherence")

    plt.tight_layout()
    plt.savefig(fname.replace(".txt","-{}.pdf".format(target)))

    if show: plt.show()

reading_data = False

with open(fname,"r") as f:
    for line in f:
        if "Coherence signal results" in line:
            reading_data = True
            continue

        if reading_data:
            if "#" in line:
                if "k_DD" in line:
                    k_DD = int(line.split()[2])
                    continue
                try: make_plot(target,k_DD,f_DD,coherence)
                except: None
                target = line.split()[2]
                f_DD = []
                coherence = []
            else:
                if nums_in(line) != 2: continue
                f_DD.append(float(line.split()[0]))
                coherence.append(float(line.split()[1]))

try: make_plot(target,k_DD,f_DD,coherence)
except: None
