#!/usr/bin/env python
import os

eigen_libs = ["/usr/include/eigen3/"]
if not all([ os.path.isdir(lib) for lib in eigen_libs ]):
    eigen_libs = ["~/.local/include/","~/.local/include/eigen3/"]

eigen_dirs = " ".join(["-I " + os.path.expanduser(lib) for lib in eigen_libs])
with open(".eigen-dirs","w") as f:
    f.write(eigen_dirs+"\n")
