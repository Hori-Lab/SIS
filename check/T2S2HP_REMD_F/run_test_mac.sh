mpiexec -n 4 ./mpirun.sh
mkdir test
mv out* err* test_REMD* test

mpiexec -n 4 ./mpirun1.sh
mkdir test1
mv out* err* test_REMD* test1

cat test1/test_REMD1_0*.rst > test1/test_REMD1.rst

mpiexec -n 4 ./mpirun2.sh
mkdir test2
mv out* err* test_REMD* test2
