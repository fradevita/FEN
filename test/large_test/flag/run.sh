#!/bin/bash
echo "    running flapping flag test case ..."

mkdir -p .fmod .fobj data
make SOURCE=main.f90 > compilation_log 2> compilation_warnings
mpirun -n 8 ./run.e > run_log
