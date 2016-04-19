#!/usr/bin/env python3
import sys, matplotlib, numpy
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) not in [2,3]:
    print('useage: %s file [show]' % sys.argv[0])
    exit(1)

scan_fname = sys.argv[1]
show = (len(sys.argv) == 3)

if 'scan' not in scan_fname:
    print('invalid input file')
    exit(1)

w_scan, coherence = np.loadtxt(scan_fname,unpack=True)
w_larmor, A_perp = np.loadtxt(scan_fname.replace('scan','larmor'),unpack=True)

f_larmor = w_larmor/(2*np.pi*1e3) # in kHz
f_scan = w_scan/(2*np.pi*1e3) # in kHz
A_perp /= max(A_perp)

plt.plot(f_scan,coherence,'k-',linewidth=2)
for i in range(len(f_larmor)):
    dot = 1-A_perp[i]
    plt.plot([f_larmor[i],f_larmor[i]],[1,dot],'b-',linewidth=2)
    plt.plot(f_larmor[i],dot,'bo',markersize=8)

plt.xlim(f_scan[0],f_scan[-1])
plt.ylim(-1,1)

plt.xlabel(r'$f_{\mathrm{scan}}$ (kHz)')
plt.ylabel('Coherence')

plt.tight_layout()
plt.savefig(scan_fname.replace('.txt','.pdf'))

if show:
    plt.show()
