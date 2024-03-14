#!/bin/bash

echo "Running 3D lid driven cavity test case ..."

make SOURCE=main.f90 > compilation_log 2> compilation_warning
 
mpirun -n 8 ./run.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
