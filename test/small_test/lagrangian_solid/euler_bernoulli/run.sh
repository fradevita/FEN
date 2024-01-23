#!/bin/bash

mkdir -p ./.fmod ./.fobj data
echo "    running Euler-Bernoulli test ... "
make SOURCE=euler_bernoulli.f90 > compilation_log 2> compilation_warnings
./run.e > log 2> error
python3 postpro.py $1
