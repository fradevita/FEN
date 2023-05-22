#!/run/bash

echo "Running test case for the Fast Direct Solver of the Poisson equation ..."

mkdir -p ./.fmod ./.fobj
make SOURCE=test.f90 > test_compilation_log 2> test_compilation_warnings
echo "    case 1 ..."
mpirun -n 1 ./run.e 1
echo "    case 2 ..."
mpirun -n 1 ./run.e 2
echo "    case 3 ..."
mpirun -n 1 ./run.e 3

make clean > /dev/null
make CPPDEFS+=-DDIM=3 SOURCE=test.f90 > test_compilation_log 2> test_compilation_warnings
echo "    3D case ..."
mpirun -n 1 ./run.e 2
python3 postpro.py $1

