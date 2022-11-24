#!/bin/bash

echo "Running isotropic turbulence flow test case ..."

make SOURCE=isotropic > compilation_log 2> compilation_warning
 
mpirun -n 8 ./code.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
