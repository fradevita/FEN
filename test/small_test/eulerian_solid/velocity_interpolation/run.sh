#!/bin/bash

echo "    running 2D eulerian velocity interpolation test case ..."

mkdir -p data .fmod .fobj

make SOURCE=main.f90 > compilation_log 2> compilation_warning
mpirun -n 1 ./run.e > log 2> error.err
python3 postpro.py $1

make clean > /dev/null
echo "    running 3D eulerian velocity interpolation test case ..."

make CPPDEFS+=-DDIM=3 SOURCE=main.f90 > compilation_log 2> compilation_warning
mpirun -n 1 ./run.e > log 2> error.err
python3 postpro3D.py $1

