#!/usr/bin/env python3

import os, sys, subprocess

if len(sys.argv) != 2:
    print("please specify a number of searches to perform")

searches = int(sys.argv[1])

with open(os.devnull, 'w') as null:
    subprocess.call(['fac'],stdout=null)

    found = 0
    for s in range(searches):
        found += subprocess.call(['./simulate','--pair_search','--seed',str(s+1)],stdout=null)

print("pairs found in",found,"of",searches,"searches")
