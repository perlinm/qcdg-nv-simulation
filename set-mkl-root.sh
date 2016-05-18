#!/usr/bin/env sh

mklroot_file=".mkl-root"

if [ -n "$(cat /etc/*-release | grep 'Arch Linux')" ]; then
  MKLROOT=`find /opt -name mkl | grep linux | sort --reverse | head -n 1`
fi

if ! [[ $MKLROOT && ${MKLROOT-x} ]]; then
  echo "MKLROOT not set"
  exit 1
fi

echo $MKLROOT > $mklroot_file
