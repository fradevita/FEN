#!/bin/bash
echo "    running RBC initialization test case ..."

mkdir -p data .fmod .fobj
make SOURCE=main.f90 > compilation_log 2> compilation_warnings
./run.e > run_log
python3 postprocessing.py

