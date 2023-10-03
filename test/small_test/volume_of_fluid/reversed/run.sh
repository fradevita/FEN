#!/bin/bash

mkdir -p ./.fmod ./.fobj

echo "    Running vof reversed advection test case ..."
mkdir -p data
make SOURCE=reversed.f90 > compilation_log 2> compilation_warning
mpirun -n 4 ./run.e > log 2> error.err
python3 postpro.py $1
echo "        Check image, test completed."
