#!/bin/bash

echo "Running IO restart test case ..."
mkdir -p .fmod .fobj data

make SOURCE=main.f90 > compilation_log 2> compilation_warnings
mpirun -n 1 ./run.e > run_log 2> run_warnings
diff data/state_0000002.raw data/state_0000003.raw > diff_out
NL=$(wc diff_out | awk '{print $1}')
if [[ $NL -gt 0 ]]
then
    echo "        An error occured in the IO restart test."
else
    echo "        IO restart test completed."
fi

