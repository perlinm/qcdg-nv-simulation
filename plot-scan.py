#!/usr/bin/python3
import sys, matplotlib, numpy
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 2:
    print('useage: %s file' % sys.argv[0])
    exit(1)

fname = sys.argv[1]

if 'larmor' in fname:
    larmor_fname = fname
    scan_fname = fname.replace('larmor','scan')
elif 'scan' in fname:
    scan_fname = fname
    larmor_fname = fname.replace('scan','larmor')
else:
    print('invalid data file')
    exit(1)

w_larmor, A_perp = np.loadtxt(larmor_fname,unpack=True)
w_scan, coherence = np.loadtxt(scan_fname,unpack=True)

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
plt.show()
