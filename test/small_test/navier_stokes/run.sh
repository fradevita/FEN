#!/bin/bash

echo 'Running Navier-Stokes solver test cases ...'

# Run the test suite
for d in ./*/ ; do (cd "$d" && sh run.sh $1); done

# Check for failed test
for d in ./*/
do
    cd "$d"
    if [ -f error.err ]; then
        WC=$(wc error.err | awk '{print $1}')
        if [ "$WC" -ne "0" ]; then
            echo "ERROR in the ${d:2:-1} test case";
        fi
    fi
    cd ..
done

echo 'Navier-Stokes test cases completed.'
