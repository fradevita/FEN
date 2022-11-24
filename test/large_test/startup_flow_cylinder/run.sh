#!/bin/bash

echo "Running startup flow around a cylinder at Re = 1000 test case ..."

make SOURCE=main > compilation_log 2> compilation_warning
 
mpirun -n 8 ./code.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
