#!/bin/bash

echo "Running Poiseuille test case with power law model ..."

make SOURCE=power_law > compilation_log 2> compilation_warning
 
mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1
