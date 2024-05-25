#!/usr/bin/env python

NUCLEOTIDES = ['A', 'U', 'G', 'C', 'D']

import sys
import numpy as np

if len(sys.argv) != 3:
    print('Usage: SCRIPT [FASTA file] [output PSF file]')
    sys.exit(2)

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
    fpsf.write(f'{i_nt:10d} R {i_nt:8d}         R{s}      C4        C4      0.000000        0.0000           0 \n')

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

fpsf.write('         0 !NPHI: dihedrals\n')
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

