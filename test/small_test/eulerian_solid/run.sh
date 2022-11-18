echo "Running 2D velocity interpolation test case ..."

make SOURCE=2D_interpolation > compilation_log 2> compilation_warning
 
mpirun -n 1 ./code.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
