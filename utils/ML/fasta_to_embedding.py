#!/usr/bin/env python

NU2INT = {'A':1,
          'G':2,
          'C':3,
          'U':4,
          'D':5,
          }

import sys
import numpy as np

if len(sys.argv) != 3:
    print('Usage: SCRIPT [FASTA file] [output embedding file]')
    sys.exit(2)

seq = []
for l in open(sys.argv[1]):
    if l.startswith(('>', ';')):
        continue
        
    for s in list(l.strip()):
        if s == '*':
            continue
        if s.upper() not in NU2INT.keys():
            print('Error: Unknown nucleotide type: ', s)
            sys.exit(2)
        seq.append(s.upper())

emb = np.array([NU2INT[s] for s in seq], dtype='<U3')

np.save(sys.argv[2], emb)
