#!/usr/bin/env python

import os.path
import sys
import glob
import numpy as np

if len(sys.argv) != 3:
    print('Usage: SCRIPT [filepath prefix of force files] [prefix for output npy files]')
    print('             e.g. ~/ML/ for ~/ML/force_xxx.out')
    sys.exit(2)

filepath = sys.argv[1]
filepath_out = sys.argv[2]

force_out_files = glob.glob(filepath + 'force_*.out')

for out_file in force_out_files:

    data = []
    with open(out_file) as f:
        for line in f:
            try:
                # Attempt to convert each value in the line to a float
                values = [float(value) for value in line.strip().split()]
                data.append(values)
            except ValueError:
                # Ignore lines that cannot be converted to float
                pass

    Naa = len(data[0]) // 3
    numpy_array = np.array(data).reshape(-1, Naa, 3)

    npy_path = filepath_out + '/' + os.path.basename(out_file)[0:-4] + '.npy'
    np.save(npy_path, numpy_array)