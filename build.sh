set -eu

mkdir -p build
cd build

hostname=`hostname`
hostnameall=`hostname -A`

if [[ ${hostnameall:(-14):(-1)} == 'archer2.ac.uk' ]]; then
    #FC=ftn cmake .. -DCMAKE_BUILD_TYPE=Release -DMT_USE_NORANGECHECK_INSTEAD_OF_ALLOWINVALIDBOZ=ON
    FC=ftn cmake .. -DCMAKE_BUILD_TYPE=Release -DMT_USE_ALLOWINVALIDBOZ=OFF -DBP_HALT_IEEE_EXCEPTIONS=OFF
elif [[ ${hostname:8:9} == 'sulis.hpc' ]]; then
    FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
elif [[ ${hostname:13:24} == 'augusta.nottingham.ac.uk' ]]; then
    FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
elif [[ ${hostname:12:5} == 'cosma' ]]; then
    FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
elif [[ ${hostname:(-22)} == 'pharm.nottingham.ac.uk' ]]; then
    FC=gfortran cmake .. -DCMAKE_BUILD_TYPE=Release -DMT_USE_NORANGECHECK=ON -DMT_USE_ALLOWINVALIDBOZ=OFF
else
    #FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug
    FC=gfortran cmake .. -DCMAKE_BUILD_TYPE=Release
fi

make -j8
#make VERBOSE=1

cd ../
ln -sf build/sis
