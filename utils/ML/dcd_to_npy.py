#!/usr/bin/env python

from dcd import DcdFile, DcdHeader

if __name__ == '__main__':

    import sys
    import numpy as np

    if len(sys.argv) != 3:
        print('Usage: SCRIPT [dcd file] [output npy file]')
        sys.exit(2)

    dcd = DcdFile(sys.argv[1])
    dcd.open_to_read()
    dcd.read_header()
    n_frame = dcd.count_frame()
    n_nt = dcd._header.nmp_real

    print("N_nt = ", n_nt)
    print("N_frame = ", n_frame)

    data = []
    while dcd.has_more_data():
        xyz = dcd.read_onestep()
        data.append(xyz)

    np.save(sys.argv[-1], np.array(data))

    dcd.close()

