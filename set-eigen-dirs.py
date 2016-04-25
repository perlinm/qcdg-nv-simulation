#!/usr/bin/env python3
import os

eigen_dirs = "-I /usr/include/eigen3/"
if not os.path.isdir(eigen_dirs.split()[-1]):
    eigen_dirs = "-I ~/.local/include/ -I ~/.local/include/eigen3/"

with open(".eigen-dirs","w") as f:
    f.write(eigen_dirs+"\n")
