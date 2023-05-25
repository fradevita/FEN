#!/bin/bash

echo "Running Eulerian tagging test case ..."

make SOURCE=tagging.f90 > compilation_log 2> compilation_warning
 
mpirun -n 1 ./run.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
