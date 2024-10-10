#!/bin/bash

mkdir -p ./.fmod ./.fobj

echo "    Running Zalesak test case ..."
make SOURCE=Zalesak.f90 > compilation_log 2> compilation_warning

echo "        running case 1 ..."
mpirun -n 1 ./run.e 1 > log 2> error.err

echo "        running case 2 ..."
mpirun -n 1 ./run.e 2 > log 2> error.err

echo "        running case 3 ..."
mpirun -n 1 ./run.e 3 > log 2> error.err

python3 postpro.py $1
echo "        Check image, test completed."
