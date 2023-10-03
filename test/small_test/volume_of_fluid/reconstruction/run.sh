#!/bin/bash

mkdir -p ./.fmod ./.fobj

echo "    Running vof reconstruction test case ..."
make SOURCE=reconstruction.f90 > compilation_log 2> compilation_warning
mpirun -n 1 ./run.e > log 2> error.err
python3 postpro.py $1
echo "        Check image, test completed."
