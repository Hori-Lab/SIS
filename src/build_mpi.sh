mkdir -p ../build_mpi && cd ../build_mpi

## MPI gfortran
#FC=mpifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_MPI=ON

hostname=`hostname`

if [[ ${hostname:8:9} == 'sulis.hpc' ]]; then
    FC=mpiifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_MPI=ON
elif [[ ${hostname:13:24} == 'augusta.nottingham.ac.uk' ]]; then
    #FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
    echo 'Edit build_mpi.sh!'
    exit 2
elif [[ ${hostname:12:5} == 'cosma' ]]; then
    FC=mpiifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_MPI=ON
else
    FC=mpifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_MPI=ON
fi

make VERBOSE=1

cd ../src
