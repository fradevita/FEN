#!/bin/bash

# Convert ASCII stl to input file for the program
grep 'vertex' $1 | awk '{print $2 " " $3 " " $4}' > $2
