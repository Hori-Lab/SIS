#!/bin/sh

rank=$OMPI_COMM_WORLD_RANK

suffix=`printf "%3.3d" $rank`

../sis input_REMD.toml 1> ./out.$suffix 2> ./err.$suffix
