echo "Running Euler-Bernoulli test case..."

python3 mesh.py

make SOURCE=Euler_Bernoulli > compilation_log 2> compilation_warning

mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1
