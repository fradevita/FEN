#!/bin/bash

echo "    Running rising bubble test case..."
mkdir -p ./.fmod ./.fobj

make SOURCE=rising_bubble.f90 > compilation_log 2> compilation_warning

echo "        Running test case 1 ..."
mpirun -n 4 ./run.e 1 > log 2> error.err
python3 postpro.py $1
echo "            Check image, test completed."

echo "        Running test case 2 ..."
mpirun -n 4 ./run.e 2 > log 2> error.err
python3 postpro_2.py $1
echo "            Check image, test completed"
