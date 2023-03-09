#!/bin/bash

echo "Running laminar cylinder at Re = 20 test case with Eulerian IBM"

mkdir -p data

make SOURCE=cylinder > compilation_log 2> compilation_warning
mpirun -n 4 ./code.e > log 2> error.err

sed -i 's/Nx = 64/Nx = 128/' cylinder.f90
make SOURCE=cylinder > compilation_log 2> compilation_warning
mpirun -n 4 ./code.e > log 2> error.err

sed -i 's/Nx = 128/Nx = 256/' cylinder.f90
make SOURCE=cylinder > compilation_log 2> compilation_warning
mpirun -n 4 ./code.e > log 2> error.err

sed -i 's/Nx = 256/Nx = 64/' cylinder.f90

python3 postpro.py $1
