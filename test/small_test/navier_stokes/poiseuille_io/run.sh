#!/bin/bash

echo "    Running Poiseuille InOut test case ..."

mkdir -p ./.fmod ./.fobj data
make SOURCE=poiseuille_io.f90 > compilation_log 2> compilation_warning
mpirun -n 4 ./run.e > log 2> error.err
python3 postpro.py $1
