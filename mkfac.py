#!/usr/bin/env python3

import sys

error_flags = '-Wall -Werror '
if len(sys.argv) > 1:
    error_flags = ''

fac_text = '''| g++ {0}-O3 -std=c++14 -I /usr/include/eigen3/ -g -flto -c -o qp-math.o qp-math.cpp
< qp-math.h
< qp-math.cpp
> qp-math.o

| g++ {0}-O3 -std=c++14 -I /usr/include/eigen3/ -g -flto -c -o nv-math.o nv-math.cpp
< nv-math.h
< nv-math.cpp
> nv-math.o

| g++ {0}-O3 -std=c++14 -lboost_program_options -lboost_system -lboost_filesystem -I /usr/include/eigen3/ -g -flto -c -o simulation.o simulation.cpp
< simulation.cpp
> simulation.o

| g++ {0}-std=c++14 -lboost_program_options -lboost_system -lboost_filesystem -flto -o simulate qp-math.o nv-math.o simulation.o
< qp-math.o
< nv-math.o
< simulation.o
> simulate
'''

with open('./.simulation.fac','w') as f:
    f.write(fac_text.format(error_flags))
