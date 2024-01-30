echo "    running MLS interpolation 2D test case ..."

mkdir -p ./.fmod ./.fobj
CPPFLAGS+=-DDIM=2 make SOURCE=test_2D.f90 > compilation_log 2> compilation_warning
mpirun -n 1 ./run.e > log 2> error.err
python3 postpro_2D.py $1

echo "    running MLS interpolation 3D test case ..."
CPPFLAGS+=-DDIM=3 make SOURCE=test_3D.f90 > compilation_log 2> compilation_warning
mpirun -n 1 ./run.e > log 2> error.err
python3 postpro_3D.py $1
