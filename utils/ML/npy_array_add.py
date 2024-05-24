#!/usr/bin/env python

import sys
import numpy as np

if len(sys.argv) == 1:
    print('Usage: SCRIPT [input npy files [npy files ...]] [output npy file]')
    sys.exit(2)

files = sys.argv[1:-1]
out_npy = sys.argv[-1]

a = np.load(files[1])

for f in files[1:]:
    a += np.load(f)

np.save(out_npy, a)
