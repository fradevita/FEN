#!/bin/bash

echo "Running capillary wave test case..."

mkdir -p data_08 data_16 data_32 data_64

make SOURCE=capillary > compilation_log 2> compilation_warning

mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1
