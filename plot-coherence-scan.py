#!/usr/bin/env python
import sys, matplotlib
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) not in [2,3]:
    print('useage: %s file [show]' % sys.argv[0])
    exit(1)

fname = sys.argv[1]
output_fname = fname.replace(".txt","") + ".pdf"

show = (len(sys.argv) == 3)

reading_larmor_data = False
reading_scan_data = False

w_larmor = []
hyperfine_perp = []
w_scan = []
coherence = []

with open(fname, "r") as f:
    for line in f:
        if line[0] == "#": continue

        if "-----" in line:
            if not reading_larmor_data:
                reading_larmor_data = True
            elif not reading_scan_data:
                reading_scan_data = True
                reading_larmor_data = False
            continue

        if reading_larmor_data:
            w_larmor.append(float(line.split()[0]))
            hyperfine_perp.append(float(line.split()[1]))
            continue
        if reading_scan_data:
            w_scan.append(float(line.split()[0]))
            coherence.append(float(line.split()[1]))

f_larmor = np.array(w_larmor)/(2*np.pi*1e3)
f_scan = np.array(w_scan)/(2*np.pi*1e3)
relative_hyperfine = np.array(hyperfine_perp)/max(hyperfine_perp)

plt.plot(f_scan,coherence,'k-',linewidth=2)
for i in range(len(f_larmor)):
    dot = 1 - relative_hyperfine[i]
    plt.plot([f_larmor[i],f_larmor[i]],[1,dot],'b-',linewidth=2)
    plt.plot(f_larmor[i],dot,'bo',markersize=8)

plt.xlim(f_scan[0],f_scan[-1])
plt.ylim(-1,1)

plt.xlabel(r'$f_{\mathrm{scan}}$ (kHz)')
plt.ylabel('Coherence')

plt.tight_layout()
plt.savefig(output_fname)

if show:
    plt.show()
