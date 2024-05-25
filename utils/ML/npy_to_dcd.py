#!/usr/bin/env python

from dcd import DcdFile, DcdHeader

def move_to_com(xyz):
    return xyz - np.mean(xyz, axis=0)

if __name__ == '__main__':

    import sys
    import numpy as np

    if len(sys.argv) != 3:
        print('Usage: SCRIPT [npy file] [output dcd file]')
        sys.exit(2)

    data = np.load(sys.argv[1])
    n_nt, _, n_frame = np.shape(data)
    print("N_nt = ", n_nt)
    print("N_frame = ", n_frame)

    dcd = DcdFile(sys.argv[-1])
    dcd.open_to_write()

    # Prepare a header
    header = DcdHeader()
    header.title = []
    header.nset = n_frame
    header.nmp_real = n_nt

    # Write the header
    dcd.set_header(header)
    dcd.write_header()

    # Write the trajectory
    for d in data.T:
        dcd.write_onestep(move_to_com(d.T))

    dcd.close()

