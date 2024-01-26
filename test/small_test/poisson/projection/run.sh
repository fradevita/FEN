#!/run/bash

echo "    Running test case for the projection operator ..."

mkdir -p ./.fmod ./.fobj
make SOURCE=projection.f90 > compilation_log 2> compilation_warnings
mpirun -n 1 ./run.e 1 > log
mpirun -n 1 ./run.e 2 >> log

NL=$(wc log | awk '{print $1}')
if [[ $NL -gt 0 ]]
then
    echo "        An error occured in the projecton test."
else
    echo "        projection test completed."
fi

