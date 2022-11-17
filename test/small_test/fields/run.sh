echo "Running test field ..."

make SOURCE=test_fields > compilation_log 2> compilation_warning

mpirun -n 8 ./code.e > log 2> error.err

NL=$(wc error.err | awk '{print $1}')
if [[ $NL -gt 0 ]]
then
    echo "An error occured in the test field."
else
    echo "Test completed."
fi
