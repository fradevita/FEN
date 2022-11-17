echo "Running Poiseuille test case ..."

make SOURCE=test_Poiseuille > compilation_log 2> compilation_warning
 
mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1
