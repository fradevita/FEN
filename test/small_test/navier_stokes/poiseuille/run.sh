#!/bin/bash

echo "    Running Poiseuille test case ..."

make SOURCE=poiseuille.f90 > compilation_log 2> compilation_warning
 
mpirun -n 1 ./run.e > log 2> error.err

python3 postpro.py $1
