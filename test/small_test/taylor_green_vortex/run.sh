echo "Running Taylor Green Vortex test case ..."

make SOURCE=test_TGV > compilation_log 2> compilation_warning
 
mpirun -n 4 ./code.e > log 2> error.err

python3 postpro.py $1
