echo "Running rising bubble test case..."

make SOURCE=rising_bubble > compilation_log 2> compilation_warning

mpirun -n 4 ./code.e > log 2> error.err

python3 postpro.py $1

echo "Check image, test completed."
