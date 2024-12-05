#!/bin/bash

echo "Running IO restart test case ..."
mkdir -p .fmod .fobj data

CPPDEFS+=-DDIM=3 make SOURCE=test_NS.f90 > compilation_log 2> compilation_warnings
mpirun -n 1 ./run.e > run_log 2> run_warnings
diff data/state_0000010.raw data/state_0000020.raw > diff_out
NL=$(wc diff_out | awk '{print $1}')
if [[ $NL -gt 0 ]]
then
    echo "        An error occured in the IO NS restart test."
else
    echo "        IO NS restart test completed."
fi

bash clean.sh > /dev/null
CPPDEFS+=-DMF make SOURCE=test_MF.f90 > compilation_log 2> compilation_warnings
mpirun -n 1 ./run.e > run_log 2> run_warnings
diff data/state_0000010.raw data/state_0000020.raw > diff_out
NL=$(wc diff_out | awk '{print $1}')
if [[ $NL -gt 0 ]]
then
    echo "        An error occured in the IO MF restart test."
else
    echo "        IO MF restart test completed."
fi

bash clean.sh > /dev/null
CPPDEFS+=-DIBM CPPDEFS+=-DFSI make SOURCE=test_FSI_E.f90 > compilation_log 2> compilation_warnings
mpirun -n 1 ./run.e > run_log 2> run_warnings
diff data/state_0000010.raw data/state_0000020.raw > diff_out
NL=$(wc diff_out | awk '{print $1}')
if [[ $NL -gt 0 ]]
then
    echo "        An error occured in the IO FSI_E restart test."
else
    echo "        IO FSI_E restart test completed."
fi
