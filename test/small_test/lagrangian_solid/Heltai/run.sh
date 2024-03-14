#!/bin/bash

mkdir -p data
echo "    running 2D cantilever test ... "
make SOURCE=main.f90 > compilation_log 2> compilation_warnings
./run.e > run_log
