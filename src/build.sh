mkdir -p ../build && cd ../build

#cmake -DCMAKE_BUILD_TYPE=Debug ..
#cmake -DCMAKE_BUILD_TYPE=Release ..

# Fortran debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-O0 -g -fbacktrace -ffpe-trap=zero,overflow,denormal,underflow -Wall -ffree-line-length-none -fbounds-check -fmax-errors=5" ..

#cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS="-heap-arrays" ..
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-O0 -g -traceback -check all -heap-arrays" ..

make -j12
#make -j12 VERBOSE=1

cd ../src
