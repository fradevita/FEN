#!/bin/bash

echo 'Running memory check test case for lagrangian solid ...'

make clean > /dev/null 
make SOURCE=mem_test > compilation_log 2> compilation_warning
python3 mesh.py
valgrind --leak-check=full --show-leak-kinds=all ./code.e
