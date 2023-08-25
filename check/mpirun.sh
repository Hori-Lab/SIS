#!/bin/sh

rank=$OMPI_COMM_WORLD_RANK
#rank=$PMI_RANK

suffix=`printf "%3.3d" $rank`

../sis.mpi input_REMD.toml 1> ./out.$suffix 2> ./err.$suffix
