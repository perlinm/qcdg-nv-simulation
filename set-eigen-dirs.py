#!/usr/bin/env python
import os

eigen_dirs = "-I /usr/include/eigen3/"
if not os.path.isdir(eigen_dirs.split()[-1]):
    eigen_libs = ["~/.local/include/","~/.local/include/eigen3/"]
    eigen_dirs = " ".join(["-I " + os.path.expanduser(lib) for lib in eigen_libs])

with open(".eigen-dirs","w") as f:
    f.write(eigen_dirs+"\n")
