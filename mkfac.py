#!/usr/bin/env python3
import sys, os, glob, re

executable = "simulate"
sim_files = sorted(glob.glob("*.cpp"))

std = "-std=c++11"
debug_info = "-g"
optimization = "-O3"
error_flags = "-Wall -Werror"
testing_mode = (len(sys.argv) > 1)

mkl_root = ".mkl-root"
eigen_libs = ".eigen-libs"

mkl_flags = ("-Wl,--no-as-needed,-rpath=`cat {0}`/lib/intel64/ -L `cat {0}`/lib/intel64/" + \
             " -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl -fopenmp" + \
             " -m64 -I `cat {0}`/include/").format(mkl_root)

lib_flags = {"eigen3" : "`cat {}`".format(eigen_libs),
             "boost/filesystem" : "-lboost_system -lboost_filesystem",
             "boost/program_options" : "-lboost_program_options"}

global_dependencies = [ mkl_root, eigen_libs ]

fac_text = ""
all_libraries = []
all_headers = []

def fac_rule(libraries, headers, out_file, in_files, link=False):
    text = "| g++ {} -flto ".format(std)
    if not link: text += "-c "
    if testing_mode:
        text += debug_info + " "
    else:
        text += "{} {} {} ".format(optimization,error_flags,mkl_flags)
    text += " ".join(libraries)
    text += " -o {} ".format(out_file)
    text += " ".join(in_files)+"\n"
    for dependency in headers + in_files + global_dependencies:
        text += "< {}\n".format(dependency)
    text += "> {}\n\n".format(out_file)
    return text

for sim_file in sim_files:
    out_file = sim_file.replace(".cpp",".o")
    libraries = []
    headers = []
    with open(sim_file,'r') as f:
        for line in f:
            if "#include" in line:
                for flag in lib_flags.keys():
                    if flag in line and lib_flags[flag] not in libraries:
                        libraries += [lib_flags[flag]]
                if re.search('"*.h"',line):
                    headers += [line.split('"')[-2]]


    fac_text += fac_rule(libraries, headers, out_file, [sim_file])
    for library in libraries:
        if library not in all_libraries:
            all_libraries += [library]
    for header in headers:
        if header not in all_headers:
            all_headers += [header]


out_files = [ sim_file.replace(".cpp",".o") for sim_file in sim_files ]
fac_text += fac_rule(all_libraries, all_headers, executable, out_files, link=True)

with open(".{}.fac".format(executable),"w") as f:
    f.write(fac_text)
