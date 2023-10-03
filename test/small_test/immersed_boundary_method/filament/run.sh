echo "    running filament in fluid test case..."

python3 mesh.py
make SOURCE=filament.f90 > compilation_log 2> compilation_warning
mpirun -n 8 ./run.e > log 2> error.err
python3 postpro.py $1
