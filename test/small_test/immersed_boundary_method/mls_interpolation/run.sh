echo "    Running MLS interpolation test case ..."

mkdir -p ./.fmod ./.fobj
make SOURCE=test.f90 > compilation_log 2> compilation_warning
mpirun -n 4 ./run.e > log 2> error.err
python3 postpro.py $1
