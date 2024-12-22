#!/usr/bin/env python

import math
import argparse

parser = argparse.ArgumentParser(
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('outfile', type=argparse.FileType('w'),
                    help='Output filename. The output format depends on the extension, either .pdb or .xyz.')

group_seq = parser.add_mutually_exclusive_group(required=True)
group_seq.add_argument('--seq', help='Sequence')
group_seq.add_argument('--fasta', type=argparse.FileType('r'), help='Sequence FASTA file')

parser.add_argument('--b', type=float, default=6.13, help='Distance between adjacent nucleotides')

args = parser.parse_args()

if args.seq is not None:
    seq = args.seq
    n_nt = len(seq)

elif args.fasta is not None:
    seq = ''
    for l in args.fasta:
        if l[0] == '>' or l[0] == '#' or l[0] == ';':
            continue
        seq += l.strip()
    n_nt = len(seq)

n_nt = len(seq)
print ('Sequence: %s' % (seq,))
print ('#nucleotides: %i' % (n_nt,))
b = args.b   #5.5  # Distance between P and P

theta = 2 * math.pi / float(n_nt)
r = b / theta

if args.outfile.name[-3:] == 'xyz':
    args.outfile.write(f'{n_nt}\n')
    args.outfile.write(f'\n')

for i in range(n_nt):
    s = seq[i-1]

    x = r * math.cos(theta * i)
    y = 0
    z = r * math.sin(theta * i)
    if args.outfile.name[-3:] == 'xyz':
        args.outfile.write('%s %8.3f %8.3f %8.3f\n' % (s,x,y,z))
    else:
        args.outfile.write('ATOM  %5i  P   R%s  A %3i    %8.3f%8.3f%8.3f\n' % (i+1, s, i+1, x,y,z))

args.outfile.close()
