#!/bin/bash
echo "    running 2D stretching test case ..."

basedir="data_"
underscore="_"
\rm -rf data_*
for arg2 in {0..2}
do
    for arg1 in {0..5}
    do
        datadir="$basedir$arg1$underscore$arg2"
        echo "        running case with " $arg1 $arg2
        make SOURCE=main.f90 > compilation_log 2> compilation_warnings
        mkdir -p data
        ./run.e $arg1 $arg2# > run_log
        mv data $datadir
    done
done

python3 postprocessing.py $1
