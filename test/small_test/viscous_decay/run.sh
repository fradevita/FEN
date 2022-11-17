echo "Running viscous decay test case..."

make SOURCE=viscous_decay > compilation_log 2> compilation_warning

mpirun -n 4 ./code.e > log 2> error.err

python3 postpro.py $1
