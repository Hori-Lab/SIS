#!/usr/bin/env python

import numpy as np
import argparse

## Parameters
Nback = 10   # Number of particles removed backward when new bead causes clash.
BL = 5.84    # Bond length
BA = 2.643   # Bond angle
EV = 12.0    # Excluded-volume distance
EV2 = EV * EV

# Arguments parser
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--xyz', default='random.xyz', help='output XYZ filename')
parser.add_argument('--nxyz', help='output XYZ filename for Kymoknot')
group_seq = parser.add_mutually_exclusive_group(required=True)
group_seq.add_argument('--fasta', type=argparse.FileType('r'), help='Sequence FASTA file')
group_seq.add_argument('--Nbeads', type=int, help='Number of beads')

args = parser.parse_args()

if args.Nbeads is not None:
    N = args.Nbeads

seq = ''
if args.fasta is not None:
    for l in args.fasta:
        if l.startswith(('>', ';')):
            continue
        if len(l.strip()) == 0:
            continue
        seq += l.strip()

    N = len(seq)

else:
    seq = 'A'*N

# Natural Extension of Reference Frame
def NeRF(A, B, C, bond, angl, dihd):

    t = np.pi - angl

    D2 = np.array([bond * np.cos(t),
                   bond * np.cos(dihd) * np.sin(t),
                   bond * np.sin(dihd) * np.sin(t)])
  
    AB = B - A
    BC = C - B
    n_bc = BC / np.linalg.norm(BC)

    n = np.cross(AB, n_bc)
    n = n / np.linalg.norm(n)

    Mx = n_bc
    My = np.cross(n, n_bc)
    #My = My / np.linalg.norm(My)   # Not necessary as the norm should be 1
    Mz = n
    M = np.array([Mx, My, Mz]).T
    # Apply transpose becasue vectors (Mx, My, Mz) have to be the columns of M.

    D = np.matmul(M, D2) + C

    return D

# Coordinates for the first three beads.
xyz = []
xyz.append(np.array([0., 0., 0.]))
xyz.append(np.array([BL, 0., 0.]))
xyz.append(np.array([BL*(1.+np.cos(np.pi-BA)), BL*np.sin(np.pi-BA), 0.]))
i = 3

# Main loop
while (i < N):
    # Generate a dihedral angle at random
    dih = -np.pi + np.random.random() * 2 * np.pi

    # Derive the coordinate by NeRF
    p = NeRF(xyz[-3], xyz[-2], xyz[-1], BL, BA, dih)

    # Check if any clash
    flg = False
    for j in range(i-2):
        v = p - xyz[j]
        d2 = v[0]**2 + v[1]**2 + v[2]**2
        if d2 < EV2:
            flg = True
            break

    # If clash, remove the last Nback beads
    if flg:
        Ndel = min(i-3, Nback)
        i -= Ndel
        del xyz[-Ndel:]
        continue

    # Add the new bead
    xyz.append(p)
    i += 1

# Output XYZ file
f = open(args.xyz, 'w')
f.write(f'{N}\n')
f.write('\n')
for i, p in enumerate(xyz):
    f.write(f'{seq[i]}  {p[0]:10.3f}  {p[1]:10.3f}  {p[2]:10.3f}\n')
f.close()

# Output NXYZ file (optional)
if args.nxyz is not None:
    f = open(args.nxyz, 'w')
    f.write(f'{N}\n')
    for p in xyz:
        f.write(f' {p[0]:10.3f}  {p[1]:10.3f}  {p[2]:10.3f}\n')
    f.close()
