#!/bin/bash

echo "    Running lid driven test case ..."

mkdir -p ./.fmod ./.fobj data
make SOURCE=lid_driven.f90 > compilation_log 2> compilation_warning 
mpirun -n 4 ./run.e > log 2> error.err
python3 postpro.py $1

echo "        Check image, test completed."
