#!/run/bash

echo "Running test case for the Fast Direct Solver of the Poisson equation ..."

mkdir -p ./.fmod ./.fobj
make SOURCE=projection.f90 > compilation_log 2> compilation_warnings
mpirun -n 1 ./run.e
python3 postpro.py $1

