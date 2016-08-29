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

def figname(targets):
    return fname.replace(".txt","-{}.pdf".format("-".join(targets)))

def line_plot(targets,k_DD,f_DD,coherence):

    plot = plt.plot(f_DD,coherence,'k-',linewidth=2)

    plt.xlim(0,f_DD[-1])
    plt.ylim(-1,1)

    plt.xlabel("$f_{}$".format(k_DD))
    plt.ylabel("Coherence")

    plt.tight_layout()
    plt.savefig(figname(targets))

    if show: plt.show()

def color_plot(targets,k_DD,f_DD,coherence):

    f_DD_boundaries = ( f_DD - (f_DD[1]-f_DD[0])/2 )
    d_phi = 1/coherence.shape[1]
    angles_over_pi = np.arange(0,1+d_phi,1/coherence.shape[1])

    plt.pcolor(angles_over_pi,f_DD_boundaries,coherence)

    plt.xlim(0,1)
    plt.ylim(0,f_DD_boundaries[-1])
    plt.clim(-1,1)

    plt.xlabel(r"$\phi_{DD}/\pi$")
    plt.ylabel("$f_{}$".format(k_DD))

    cbar = plt.colorbar()
    cbar.set_label("Coherence")

    plt.tight_layout()
    plt.savefig(figname(targets))

    if show: plt.show()

def make_plot(targets,k_DD,f_DD,coherence):
    f_DD = np.array(f_DD)
    coherence = np.array(coherence)

    if coherence.shape[0] == 0: return None
    if coherence.shape[1] == 1:
        line_plot(targets,k_DD,f_DD,coherence)
    else:
        color_plot(targets,k_DD,f_DD,coherence)

reading_data = False

with open(fname,"r") as f:
    for line in f:
        if "Coherence signal results" in line:
            reading_data = True
            continue

        if reading_data:
            if "#" in line:
                if "k_DD" in line:
                    k_DD = int(line.split()[-1])
                    continue
                if "targets" in line:
                    try: make_plot(targets,k_DD,f_DD,coherence)
                    except: None
                    targets = line.split()[2:]
                    f_DD = []
                    coherence = []
            else:
                f_DD.append(float(line.split()[0]))
                coherence.append([ float(n) for n in line.split()[1:]])

make_plot(targets,k_DD,f_DD,coherence)
