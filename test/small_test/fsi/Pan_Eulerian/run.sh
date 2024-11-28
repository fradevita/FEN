echo "    running Pan test case ..."

mkdir -p ./.fmod ./.fobj data
make SOURCE=Pan.f90 > compilation_log 2> compilation_warning
python3 mesh.py
for i in {1..3}
do
    echo "        running case ${i}"
    mpirun -n 4 ./run.e $i > log 2> error.err
done
python3 postpro.py $1
