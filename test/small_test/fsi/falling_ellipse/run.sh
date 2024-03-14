echo "    running falling ellipse test case..."

mkdir -p data ./.fmod ./.fobj
python3 mesh.py
make SOURCE=ellipse.f90 > compilation_log 2> compilation_warning
mpirun -n 4 ./run.e > log 2> error.err
python3 postpro.py $1
