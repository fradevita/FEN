#!/bin/bash

echo 'Running norm from distance test case ...'

make SOURCE=main > compilation_log 2> compilation_warning
./code.e > log 2> error.err
python3 postpro.py $1 2> plot.err

echo 'Done'
