#!/bin/bash

mkdir -p ./.fmod ./.fobj data
echo "    running energy test ... "
make SOURCE=energy.f90 > energy_compilation_log 2> energy_compilation_warnings
python3 mesh.py
./run.e > energy_log 2> error.err
python3 postpro.py $1
