echo "    running Pan test case with lagrangian solid ..."

make SOURCE=Pan.f90 > compilation_log 2> compilation_warning
python3 mesh.py
mpirun -n 1 ./run.e > log 2> error.err
python3 postpro.py $1
