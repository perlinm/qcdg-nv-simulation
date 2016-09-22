#!/usr/bin/env python
import sys, matplotlib
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) not in [2,3,4]:
    print("useage: %s file [show] [azimuths]" % sys.argv[0])
    exit(1)

fname = sys.argv[1]

show = "show" in sys.argv
show_azimuths = "azimuths" in sys.argv

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
    suffix = "-{}".format("-".join(targets))
    if show_azimuths: suffix += "-azimuths"
    return fname.replace(".txt",suffix+".png")

def line_plot(targets,k_DD,f_DD,coherence):

    plot = plt.plot(f_DD,coherence,"k-",linewidth=2)

    plt.xlim(0,f_DD[-1])
    plt.ylim(-1,1)

    plt.xlabel("$f_{}$".format(k_DD))
    plt.ylabel("Coherence")

    plt.tight_layout()
    plt.savefig(figname(targets))

    if show: plt.show()

def mod(x,m=1):
    while x < 0: x += m
    while x >= m: x -= m
    return x

def color_plot(targets,k_DD,f_DD,coherence,azimuths):

    f_DD_boundaries = ( f_DD - (f_DD[1]-f_DD[0])/2 )
    d_phi = 1/coherence.shape[1]
    angles_over_pi = np.arange(0,1+d_phi,1/coherence.shape[1])

    plt.pcolor(angles_over_pi,f_DD_boundaries,coherence)
    if show_azimuths:
        for azimuth in azimuths:
            pi_phi = mod(1 - azimuth)
            plt.axvline(pi_phi,color="k",linewidth=2)

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

def make_plot(targets,k_DD,f_DD,coherence,azimuths):
    f_DD = np.array(f_DD)
    coherence = np.array(coherence)

    if coherence.shape[0] == 0: return None

    plt.figure()
    if coherence.shape[1] == 1:
        line_plot(targets,k_DD,f_DD,coherence)
    else:
        color_plot(targets,k_DD,f_DD,coherence,azimuths)

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
                if "azimuths" in line:
                    azimuths = [ float(a) for a in line.split()[2:]]
                    continue
                if "targets" in line:
                    try: make_plot(targets,k_DD,f_DD,coherence,azimuths)
                    except: None
                    targets = line.split()[2:]
                    f_DD = []
                    coherence = []
            else:
                f_DD.append(float(line.split()[0]))
                coherence.append([ float(n) for n in line.split()[1:]])

make_plot(targets,k_DD,f_DD,coherence,azimuths)
