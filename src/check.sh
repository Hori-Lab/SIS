cd ../check
#../sis input_check_force.toml
../sis input_T2S2HP.toml 1> ../src/out 2>&1
#../sis input_debug_ele.toml
cd ../src

cat out
