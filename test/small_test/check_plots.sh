#!/bin/bash

# check for failed test
for d in ./*/
do
    cd "$d"
    eom *.png
    cd ..
done
