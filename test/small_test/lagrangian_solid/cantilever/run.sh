#!/bin/bash

mkdir -p ./.fmod ./.fobj data
echo "    running 1D cantilever test ... "
make SOURCE=cantilever.f90 > compilation_log 2> compilation_warnings
python3 mesh.py
./run.e > log 2> error.err
python3 postpro.py $1
