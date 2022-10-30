mkdir -p ../build_mpi && cd ../build_mpi

# MPI gfortran
FC=mpifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_MPI=ON

make VERBOSE=1

cd ../src
