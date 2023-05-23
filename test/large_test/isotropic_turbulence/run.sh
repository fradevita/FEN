#!/bin/bash

echo "Running isotropic turbulence flow test case ..."

make SOURCE=isotropic.f90 > compilation_log 2> compilation_warning
 
mpirun -n 8 ./run.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
