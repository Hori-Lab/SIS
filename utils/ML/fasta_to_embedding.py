#!/usr/bin/env python

NU2INT = {'A':1,
          'G':2,
          'C':3,
          'U':4,
          'D':5,
          }

import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(
         description='Convert FASTA format file to embedding for TorchMD',
         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('fastafile', help='input FASTA file')
parser.add_argument('outfile', nargs='?', default='embedding.npy', help='output npy filename')

parser.add_argument('--noDummy', action="store_true", help='For not using dummy(D) beads at the ends.')

args = parser.parse_args()

seq = []
for l in open(args.fastafile):
    if l.startswith(('>', ';')):
        continue

    for s in list(l.strip()):
        if s == '*':
            continue
        if s.upper() not in NU2INT.keys():
            print('Error: Unknown nucleotide type: ', s)
            sys.exit(2)
        seq.append(s.upper())

# Change the terminal to dummy beads
if args.noDummy:
    pass
else:
    seq[0] = 'D'
    seq[-1] = 'D'

emb = np.array([NU2INT[s] for s in seq], dtype='<U3')

np.save(args.outfile, emb)
