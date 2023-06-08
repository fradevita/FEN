#!/bin/bash

mkdir -p ./.fmod ./.fobj data
echo "    running cantilever test ... "
make SOURCE=cantilever.f90 > compilation_log 2> compilation_warnings
python3 mesh.py
./run.e > log 2> error
python3 postpro.py $1
