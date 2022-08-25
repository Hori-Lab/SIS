set -eu

mkdir -p build_mpi
cd build_mpi

hostname=`hostname`

if [[ ${hostname:8:9} == 'sulis.hpc' ]]; then
    #FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
    echo 'Edit build_mpi.sh!'
    exit 2
elif [[ ${hostname:13:24} == 'augusta.nottingham.ac.uk' ]]; then
    #FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
    echo 'Edit build_mpi.sh!'
    exit 2
else
    FC=mpifort cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=ON
fi

make -j8

cd ../
ln -sf build_mpi/sis sis.mpi
