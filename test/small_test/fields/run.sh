#!/bin/bash

echo "Running field test cases ..."

mkdir -p ./.fmod ./.fobj

echo "    runinng memory test ..."
make SOURCE=memory.f90 > memor_compilation_log 2> memory_compilation_warnings
valgrind --leak-check=full --show-leak-kinds=all -s ./run.e 2> memory_log
grep 'HEAP' memory_log
grep 'in use at exit' memory_log
grep 'ERROR' memory_log


