mkdir -p ../build && cd ../build

#cmake ..
#cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS="-heap-arrays" ..
cmake -DCMAKE_BUILD_TYPE=Debug ..
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-O0 -g -traceback -check all -heap-arrays" ..

make -j12
#make -j12 VERBOSE=1

cd ../src
