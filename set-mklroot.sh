#!/usr/bin/env sh

if [ -n "$(cat /etc/*-release | grep 'Arch Linux')" ]; then
  MKLROOT=/opt/intel/compilers_and_libraries_2016.2.181/linux/mkl
fi

if ! [[ $MKLROOT && ${MKLROOT-x} ]]; then
  echo "MKLROOT not set"
  exit 1
fi

echo $MKLROOT > .mklroot
