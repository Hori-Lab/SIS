#!/usr/bin/env python

NUCLEOTIDES = ['A', 'U', 'G', 'C', 'D']
MASS = {'A': 328.212,
        'U': 305.164,
        'G': 344.212,
        'C': 304.182,
        'D': 300.000,
       }
CHARGE = {'A': -1.0,
          'U': -1.0,
          'G': -1.0,
          'C': -1.0,
          'D':  0.0,
         }

import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(
         description='Convert FASTA format file to PSF for TorchMD',
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('fastafile', help='input FASTA file')
parser.add_argument('outfile', nargs='?', default='output.psf', help='output PSF filename')

parser.add_argument('--noDummy', action="store_true", help='For not using dummy(D) beads at the ends.')

args = parser.parse_args()

seq = []
for l in open(sys.argv[1]):
    if l.startswith(('>', ';')):
        continue

    for s in list(l.strip()):
        if s == '*':
            continue
        if s.upper() not in NUCLEOTIDES:
            print('Error: Unknown nucleotide type: ', s)
            sys.exit(2)
        seq.append(s.upper())

# Change the terminal to dummy beads
if args.noDummy:
    pass
else:
    seq[0] = 'D'
    seq[-1] = 'D'

#print(seq)
n_nt = len(seq)

fpsf = open(sys.argv[2], 'w')
fpsf.write('PSF NAMD\n')
fpsf.write('\n')
fpsf.write('         0 !NTITLE\n')
fpsf.write('\n')

fpsf.write(f'{n_nt:10d} !NATOM\n')
for i, s in enumerate(seq):
    i_nt = i + 1
    fpsf.write(f'{i_nt:10d} R {i_nt:8d}         R{s}      C4        R{s} {CHARGE[s]:13.6f} {MASS[s]:13.4f}           0 \n')

fpsf.write('\n')
fpsf.write('\n')
fpsf.write(f'{n_nt-1:10d} !NBOND: bonds\n')
for i in range(n_nt-1):
    fpsf.write(f' {i+1:9d} {i+2:9d}')
    if i % 4 == 3:
        fpsf.write('\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write('         0 !NTHETA: angles\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write(f'{n_nt-3:10d} !NPHI: dihedrals\n')
for i in range(n_nt-3):
    fpsf.write(f' {i+1:9d} {i+2:9d} {i+3:9d} {i+4:9d}')
    if i % 2 == 1:
        fpsf.write('\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write('         0 !NIMPHI: impropers\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write('         0 !NDON: donors\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write('         0 !NACC: acceptors\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write('         0 !NNB\n')

fpsf.write('\n')
fpsf.write('\n')

fpsf.write('         0          0 !NGRP\n')

fpsf.close()

