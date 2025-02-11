#!/bin/bash

echo "Running shear droplet test case..."

# Case 1
mkdir -p data
make CPPDEFS+='-DCASE=1' SOURCE=shear_drop.f90 > compilation_log 2> compilation_warning
echo "    Case 1 ..."
mpirun -n 8 ./run.e > log 2> error.err
mv data data_1
make clean > /dev/null

# Case 2
mkdir data
make CPPDEFS+='-DCASE=2' SOURCE=shear_drop.f90 > compilation_log 2> compilation_warning
echo "    Case 2 ..."
mpirun -n 8 ./run.e > log 2> error.err
mv data data_2
make clean > /dev/null

# Case 3
mkdir data
make CPPDEFS+='-DCASE=3' SOURCE=shear_drop.f90 > compilation_log 2> compilation_warning
echo "    Case 3 ..."
mpirun -n 8 ./run.e > log 2> error.err
mv data data_3
make clean > /dev/null

python3 deformation.py
python3 postpro.py