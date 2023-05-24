#!/bin/bash

for i in {5..8}
do
    python3 mesh.py $i
    mpirun -n 4 ./run.e $i
    #python3 postpro.py
done