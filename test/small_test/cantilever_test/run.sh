#!/bin/bash

echo 'Running cantilever test case ...'

make clean > /dev/null 
make SOURCE=cantilever > compilation_log 2> compilation_warning
python3 mesh.py
./code.e > log
python3 postpro.py
