#!/bin/bash

echo '    Runinng MLS forces test case ...'

mkdir -p ./.fobj ./.fmod
make SOURCE=forces.f90 > compilation_log 2> compilation_warnings
for i in {5..8}
do
    python3 mesh.py $i
    mpirun -n 4 ./run.e $i >> error.dat
done
python3 postpro.py
