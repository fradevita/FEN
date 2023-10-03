#!/bin/bash

# Run the test suite
for d in ./*/ ; do (cd "$d" && sh run.sh $1); done

echo 'TEST SUITE COMPLETED.'
