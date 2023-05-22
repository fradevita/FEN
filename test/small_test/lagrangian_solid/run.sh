#!/bin/bash

echo "Running lagrangian solid test cases ..."

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
make clean > /dev/null

make SOURCE=energy.f90 > energy_compilation_log 2> energy_compilation_warnings
echo "    running energy test ... "
./run.e > energy_log 2> energy_error
python3 energy_postpro.py $1

make clean > /dev/null

make SOURCE=euler_bernoulli.f90 > euler_bernoulli_compilation_log 2> euler_bernoulli_compilation_warnings
echo "    running Euler-Bernoulli test ... "
./run.e > euler_bernoulli_log 2> euler_bernoulli_error
python3 euler_bernoulli_postpro.py $1

make clean > /dev/null

make SOURCE=cantilever.f90 > cantilever_compilation_log 2> cantilever_compilation_warnings
python3 cantilever_mesh.py
echo "    running cantilevecantileverr test ... "
./run.e > euler_bernoulli_log 2> euler_bernoulli_error
python3 cantilver_postpro.py $1
