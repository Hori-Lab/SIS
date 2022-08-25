cd ../check
#../sis input_check_force.toml
../sis input_T2S2HP.toml 1> ../src/out 2>&1
#../sis input_debug_ele.toml

mpirun -n 4 ./mpirun.sh

cd ../src

cat out
