#!/bin/bash

echo "Running Rayleigh-Taylor instability test case..."

mkdir -p data

make SOURCE=Rayleigh_Taylor > compilation_log 2> compilation_warning

mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
