#!/bin/bash

echo "    Running Taylor Green Vortex test case ..."

make SOURCE=taylor_green_vortex.f90 > compilation_log 2> compilation_warning
 
mpirun -n 4 ./run.e > log 2> error.err

python3 postpro.py $1
