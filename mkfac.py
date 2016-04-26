#!/usr/bin/env python
import sys, os, glob, re

testing_mode = False
hide_warnings = False
for flag in sys.argv[1:]:
    if flag.lower() == "test":
        testing_mode = True
    elif flag.lower() == "whide":
        hide_warnings = True
    else:
        print('useage: {} [test] [whide]'.format(sys.argv[0]))
        exit(1)

executable = "simulate"
sim_files = sorted(glob.glob("*.cpp"))
ignore_dirs = [ "~/.ccache/" ]

std = "-std=c++11"
optimization = "-O3"
debug_info = "-g"
warning_flags = "-Wall -Werror"

eigen_dirs = ".eigen-dirs"
mkl_root = ".mkl-root"
mkl_flags = ("-Wl,--no-as-needed,-rpath=$(cat {0})/lib/intel64/ -L $(cat {0})/lib/intel64/" + \
             " -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl -fopenmp" + \
             " -m64 -I $(cat {0})/include/").format(mkl_root)

lib_flags = {"eigen3" : "$(cat {}) ".format(eigen_dirs) + mkl_flags,
             "boost/filesystem" : "-lboost_system -lboost_filesystem",
             "boost/program_options" : "-lboost_program_options"}

global_dependencies = [ mkl_root, eigen_dirs ]

fac_text = ""
all_libraries = []
all_headers = []

def fac_rule(libraries, headers, out_file, in_files, link=False):
    text = "| g++ {} {} -flto ".format(std, debug_info if testing_mode else optimization)
    if not hide_warnings: text += warning_flags + " "
    text += " ".join(libraries) + " "
    if not link: text += "-c "
    text += "-o {} ".format(out_file)
    text += " ".join(in_files)+"\n"
    for dependency in headers + in_files + global_dependencies:
        text += "< {}\n".format(dependency)
    for ignore_dir in ignore_dirs:
        text += "C {}\n".format(ignore_dir)
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
