#!/run/bash

echo "Running test case for the Fast Direct Solver of the Poisson equation ..."

mkdir -p ./.fmod ./.fobj
make SOURCE=convergence_rate.f90 > compilation_log 2> compilation_warnings
echo "    2D case 1 ..."
mpirun -n 1 ./run.e 1
echo "    2D case 2 ..."
mpirun -n 1 ./run.e 2
echo "    2D case 3 ..."
mpirun -n 1 ./run.e 3

make clean > /dev/null
make CPPDEFS+=-DDIM=3 SOURCE=convergence_rate.f90 > compilation_log 2> compilation_warnings
echo "    3D case 1 ..."
mpirun -n 1 ./run.e 1
echo "    3D case 2 ..."
mpirun -n 1 ./run.e 2
python3 postpro.py $1

