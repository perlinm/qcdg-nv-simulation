#!/usr/bin/env python3

import os, sys, subprocess

if len(sys.argv) != 2:
    print("please specify a number of searches to perform")

searches = int(sys.argv[1])

devnull = open(os.devnull, 'w')
subprocess.call(['fac'],stdout=devnull)
devnull.close()

found = 0
for s in range(searches):
    found += subprocess.call(['./simulate','--pair_search','--seed',str(s+1)])

print("pairs found in",found,"of",searches,"searches")
