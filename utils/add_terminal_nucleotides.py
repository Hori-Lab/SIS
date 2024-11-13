#!/usr/bin/env python

import sys
import math
import numpy as np
import argparse
from itertools import product
from nerf import NeRF

## Parameters
BL = 5.84    # Bond length
BA = 2.643   # Bond angle
DIH = 0.267
EV = 12.0    # Exclusion distance to check clash
EV2 = EV * EV

# Try angle up to +/- 20 degree
BA_trials = [BA,]
for degree in range(1, 21):
    for parity in (1, -1):
        BA_trials.append(BA + degree/180.0*math.pi * parity)

# Try dihedral adding +/- 10 degree
DIH_trials = [DIH,]
for degree in range(10, 180, 10):
    for parity in (1, -1):
        DIH_trials.append(DIH + degree/180.0*math.pi * parity)

# Arguments parser
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('input', help='input XYZ file')
parser.add_argument('output', help='output XYZ file')

parser.add_argument('--end5', default=None, help="Nucleotide letter ('A', 'U', 'G' etc) to be added at the 5' end.")
parser.add_argument('--end3', default=None, help="Nucleotide letter ('A', 'U', 'G' etc) to be added at the 3' end.")

args = parser.parse_args()

if args.end5 is None and args.end3 is None:
    parser.error('At least one of the --end5 or --end3 options is required.')

# Fead input xyz file
seq = ''
xyz = []
for il, l in enumerate(open(args.input)):
    if il == 0:
        #N = int(l)
        continue
    elif il == 1:
        continue
    else:
        lsp = l.split()
        s = lsp[0]
        x = float(lsp[1])
        y = float(lsp[2])
        z = float(lsp[3])
        seq += s
        xyz.append(np.array([x,y,z]))

N = len(xyz)

# Add a nucleotide to the 5' end
flg_added = False
if args.end5 is not None:

    for d, a in product(DIH_trials, BA_trials):

        p = NeRF(xyz[2], xyz[1], xyz[0], BL, a, d)

        # Check if any clash
        flg = False
        for j in range(2, N):
            v = p - xyz[j]
            d2 = v[0]**2 + v[1]**2 + v[2]**2
            if d2 < EV2:
                flg = True
                break

        if not flg:
            xyz.insert(0, p)
            seq = args.end5 + seq
            flg_added = True
            break

if not flg_added:
    print("Error: Failed to add the nucleotide to the 5' end.")
    sys.exit(2)

# Add a nucleotide to the 3' end
flg_added = False
if args.end3 is not None:

    for d, a in product(DIH_trials, BA_trials):

        p = NeRF(xyz[-3], xyz[-2], xyz[-1], BL, a, d)

        # Check if any clash
        flg = False
        for j in range(0, N-2):
            v = p - xyz[j]
            d2 = v[0]**2 + v[1]**2 + v[2]**2
            if d2 < EV2:
                flg = True
                break

        if not flg:
            xyz.append(p)
            seq = seq + args.end3
            flg_added = True
            break

if not flg_added:
    print("Error: Failed to add the nucleotide to the 3' end.")
    sys.exit(2)


# Output XYZ file
N = len(xyz)
f = open(args.output, 'w')
f.write(f'{N}\n')
f.write('\n')
for i, p in enumerate(xyz):
    f.write(f'{seq[i]}  {p[0]:10.3f}  {p[1]:10.3f}  {p[2]:10.3f}\n')
f.close()
