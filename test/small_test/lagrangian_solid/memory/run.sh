#!/bin/bash

mkdir -p ./.fmod ./.fobj

echo "    runinng memory test ..."
make SOURCE=memory.f90 > compilation_log 2> compilation_warnings
valgrind --leak-check=full --show-leak-kinds=all -s ./run.e 2> log
NML=$(grep 'in use at exit' log | awk '{print $6}')
NME=$(grep 'ERROR' log | awk '{print $4}')
if [[ $NML -gt 0 ]] || [[ $NME -gt 0 ]];
then
    echo "        An error occured in the memory test."
    echo "            in use at exit: ${NML}"
    echo "                     ERROR: ${NME}"
else
    echo "        memory test completed."
fi
