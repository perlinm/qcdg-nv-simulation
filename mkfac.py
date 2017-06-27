#!/usr/bin/env python3
import sys, os, glob, re, subprocess
from collections import OrderedDict

testing_mode = False
hide_warnings = False
build = True
for arg in sys.argv[1:]:
    if arg.lower() == "test":
        testing_mode = True
    elif arg.lower() == "whide":
        hide_warnings = True
    elif arg.lower() == "no-build":
        build = False
    else:
        print("useage: {} [test] [whide] [no-build]".format(sys.argv[0]))
        exit(1)

executable = "simulate.exe"
eigen_dirs = ".eigen-dirs"
boost_dir = ".boost-dir"
mkl_root = ".mkl-root"

ignore_dirs = [ "~/.ccache/" ]

language_standard_flag = "-std=c++11"
warning_flags = "-Wall -Werror"
link_time_optimization_flag = "-flto"
common_flags = [ language_standard_flag,
                 warning_flags,
                 link_time_optimization_flag ]

debug_flag = "-g"
optimization_flag = "-O3"
ignored_warning_flags = [ "-Wno-unused-variable",
                          "-Wno-unused-but-set-variable",
                          "-Wno-unused-local-typedefs" ]

mkl_flags = ("-Wl,--no-as-needed,-rpath=$(cat {0})/lib/intel64/" + \
             " -L $(cat {0})/lib/intel64/ -lmkl_intel_lp64 -lmkl_core" + \
             " -lmkl_gnu_thread -lpthread -lm -ldl -fopenmp -m64" + \
             " -I $(cat {0})/include/").format(mkl_root)

lib_flags = OrderedDict()
lib_flags["eigen3"] = ["$(cat {})".format(eigen_dirs), eigen_dirs ]
lib_flags["USE_MKL"] = ["{}".format(mkl_flags), mkl_root]
lib_flags["boost"] = ["-lboost_system"] + \
                     ([boost_dir] if os.path.isfile(boost_dir) else [])
lib_flags["boost/filesystem"] = ["-lboost_filesystem"]
lib_flags["boost/program_options"] = ["-lboost_program_options"]
lib_flags["gsl"] = ["-lgsl"]

fac_text = ""
global_libraries = []
global_dependencies = []
cpp_files = sorted(glob.glob("*.cpp"))

def fac_rule(libraries, file_dependencies, input_files, output_file, link = False):
    cmd_parts = ["| g++"]
    cmd_parts += common_flags
    cmd_parts += [ debug_flag if testing_mode else optimization_flag ]
    if hide_warnings: cmd_parts += ignored_warning_flags
    if not link: cmd_parts += ["-c"]
    cmd_parts += ["-o {}".format(output_file)]
    cmd_parts += [" ".join(input_files)]
    cmd_parts += [" ".join(libraries)]
    cmd_parts = list(filter(None,cmd_parts))

    rule_text = " ".join(cmd_parts) + "\n"

    for dependency in file_dependencies + input_files:
        rule_text += "< {}\n".format(dependency)
    for ignore_dir in ignore_dirs:
        rule_text += "C {}\n".format(ignore_dir)
    rule_text += "> {}\n\n".format(output_file)
    return rule_text

for cpp_file in cpp_files:
    output_file = cpp_file.replace(".cpp",".o")
    libraries = []
    dependencies = []
    with open(cpp_file,'r') as f:
        for line in f:
            if "#include" in line or "#define" in line:
                for tag in lib_flags.keys():
                    if tag in line and lib_flags[tag][0] not in libraries:
                        libraries += [lib_flags[tag][0]]
                        if len(lib_flags[tag]) > 1:
                             dependencies += [lib_flags[tag][1]]
                if re.search('"*\.h"',line):
                    dependencies += [line.split('"')[-2]]

    fac_text += fac_rule(libraries, dependencies, [cpp_file], output_file)
    for library in libraries:
        if library not in global_libraries:
            global_libraries += [library]
    for dependency in dependencies:
        if dependency not in global_dependencies:
            global_dependencies += [dependency]


compiled_binaries = [ cpp_file.replace(".cpp",".o") for cpp_file in cpp_files ]
fac_text += fac_rule(global_libraries, global_dependencies,
                     compiled_binaries, executable, link = True)

fac_text += "| etags *.cpp *.h\n< {}\n> TAGS\n".format(executable)

with open(".{}".format(executable.replace(".exe",".fac")),"w") as f:
    f.write(fac_text)

if build: exit(subprocess.call(["fac"]))
