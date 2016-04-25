#!/usr/bin/env python3
import os

eigen_lib = "-I /usr/include/eigen3/"
if not os.path.isdir(eigen_lib.split()[-1]):
    eigen_lib = "-I ~/.local/include/ -I ~/.local/include/eigen3/"

with open(".eigen-libs","w") as f:
    f.write(eigen_lib+"\n")
