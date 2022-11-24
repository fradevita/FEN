#!/run/bash

echo "Running test case for the Fast Direct Solver of the Poisson equation ..."

make SOURCE=test_Poisson_FDS_2D > compilation_log 2> compilation_warning
 
mpirun -n 1 ./code.e > log 2> error.err

make clean > /dev/null

make CNNFLAGS+=-DDIM=3 SOURCE=test_Poisson_FDS_3D > compilation_log 2> compilation_warning

mpirun -n 1 ./code.e > log 2> error.er

python3 postpro.py $1
