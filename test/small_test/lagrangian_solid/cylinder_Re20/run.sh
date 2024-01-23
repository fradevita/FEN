echo "    running laminar cylinder at Re = 20 test case with Lagrangian IBM"

mkdir -p data ./.fmod ./.fobj
make SOURCE=cylinder.f90 > compilation_log 2> compilation_warning
python3 mesh.py 64
mpirun -n 4 ./run.e 64 > log 2> error.err
python3 mesh.py 128
mpirun -n 4 ./run.e 128 > log 2> error.err
python3 mesh.py 256
mpirun -n 4 ./run.e 256 > log 2> error.err
python3 postpro.py $1
