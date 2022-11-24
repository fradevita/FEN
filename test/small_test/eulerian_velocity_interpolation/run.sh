#!/bin/bash

echo "Running Eulerian velocity interpolation test case ..."

make SOURCE=main > compilation_log 2> compilation_warning
 
mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1

sh clean.sh

make CNNFLAGS+=-DDIM=3 SOURCE=main > compilation_log 2> compilation_warning

mpirun -n 1 ./code.e > log 2> error.err

python3 postpro3D.py $1

echo "Check image, test completed."
