echo "Running Pan test case ..."

make SOURCE=Pan > compilation_log 2> compilation_warning
python3 mesh.py
mpirun -n 4 ./code.e > log 2> error.err
python3 postpro.py $1
