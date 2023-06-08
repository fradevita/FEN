#!/bin/bash

mkdir -p ./.fmod ./.fobj

echo "    runinng memory test ..."
make SOURCE=memory.f90 > memory_compilation_log 2> memory_compilation_warnings
python3 mesh.py
valgrind --leak-check=full --show-leak-kinds=all -s ./run.e 2> memory_log
NML=$(grep 'in use at exit' memory_log | awk '{print $6}')
NME=$(grep 'ERROR' memory_log | awk '{print $4}')
if [[ $NML -gt 0 ]] || [[ $NME -gt 0 ]];
then
    echo "        An error occured in the memory test."
    echo "            in use at exit: ${NML}"
    echo "                     ERROR: ${NME}"
else
    echo "        memory test completed."
fi
