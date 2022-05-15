mkdir -p build
cd build

hostname=`hostname`

if [[ ${hostname:8:9} == 'sulis.hpc' ]]; then
    FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
else
    cmake .. -DCMAKE_BUILD_TYPE=Release
fi

make -j8

cd ../
ln -sf build/sis
