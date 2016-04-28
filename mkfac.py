#!/usr/bin/env python
import sys, os, glob, re, subprocess

testing_mode = False
hide_warnings = False
for arg in sys.argv[1:]:
    if arg.lower() == "test":
        testing_mode = True
    elif arg.lower() == "whide":
        hide_warnings = True
    else:
        print('useage: {} [test] [whide]'.format(sys.argv[0]))
        exit(1)

executable = "simulate.exe"
ignore_dirs = [ "~/.ccache/" ]

std_flag = "-std=c++11"
optimize_flag = "-O3"
debug_flag = "-g"
warning_flags = "-Wall -Werror"

eigen_dirs = ".eigen-dirs"
mkl_root = ".mkl-root"
mkl_flags = ("-Wl,--no-as-needed,-rpath=$(cat {0})/lib/intel64/" + \
             " -L $(cat {0})/lib/intel64/ -lmkl_intel_lp64 -lmkl_core" + \
             " -lmkl_gnu_thread -lpthread -lm -ldl -fopenmp -m64" + \
             " -I $(cat {0})/include/").format(mkl_root)

lib_flags = {"eigen3" : "$(cat {}) ".format(eigen_dirs) + mkl_flags,
             "boost/filesystem" : "-lboost_system -lboost_filesystem",
             "boost/program_options" : "-lboost_program_options"}

global_dependencies = [ mkl_root, eigen_dirs ]

fac_text = ""
used_libraries = []
used_headers = []
sim_files = sorted(glob.glob("*.cpp"))

def fac_rule(libraries, headers, out_file, in_files, link=False):
    rule_text = "| g++ {} {} -flto ".format(std_flag, (debug_flag if testing_mode
                                                       else optimize_flag))
    if not hide_warnings: rule_text += warning_flags + " "
    rule_text += " ".join(libraries) + " "
    if not link: rule_text += "-c "
    rule_text += "-o {} ".format(out_file)
    rule_text += " ".join(in_files)+"\n"
    for dependency in headers + in_files + global_dependencies:
        rule_text += "< {}\n".format(dependency)
    for ignore_dir in ignore_dirs:
        rule_text += "C {}\n".format(ignore_dir)
    rule_text += "> {}\n\n".format(out_file)
    return rule_text

for sim_file in sim_files:
    out_file = sim_file.replace(".cpp",".o")
    libraries = []
    headers = []
    with open(sim_file,'r') as f:
        for line in f:
            if "#include" in line:
                for tag in lib_flags.keys():
                    if tag in line and lib_flags[tag] not in libraries:
                        libraries += [lib_flags[tag]]
                if re.search('"*.h"',line):
                    headers += [line.split('"')[-2]]


    fac_text += fac_rule(libraries, headers, out_file, [sim_file])
    for library in libraries:
        if library not in used_libraries:
            used_libraries += [library]
    for header in headers:
        if header not in used_headers:
            used_headers += [header]


out_files = [ sim_file.replace(".cpp",".o") for sim_file in sim_files ]
fac_text += fac_rule(used_libraries, used_headers, executable, out_files, link=True)
fac_text += "| etags *.cpp *.h\n< {}\n> TAGS\n".format(executable)

with open(".{}".format(executable.replace(".exe",".fac")),"w") as f:
    f.write(fac_text)

subprocess.call(["fac"])
