#!/bin/bash
echo "    running settling sphere test case ..."

basedir="data_"
underscore="_"
\rm -rf data_*
make SOURCE=main.f90 > compilation_log 2> compilation_warnings
for arg in {1..4}
do
    echo "        running case " $arg
    mkdir -p data
    mpirun -n 8 ./run.e $arg > run_log
    mv run_log data/
    datadir="$basedir$arg"
    mv data $datadir
done
