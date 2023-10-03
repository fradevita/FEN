echo "    running free filament test case..."

python3 mesh.py

make SOURCE=filament.f90 > compilation_log 2> compilation_warning

mpirun -n 1 ./run.e > log 2> error.err

python3 postpro.py $1
