#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print('Usage: SCRIPT [input PDB] [output CG XYZ/PDB]')
    print('       The output file is xyz format by default, but PDB if the filename ends with .pdb.')
    sys.exit(2)

outtype = 'xyz'
if sys.argv[2][-3:] == 'pdb':
    outtype = 'pdb'
    
output = open(sys.argv[2], 'w')

unified_base = {'A': 'A', 'U': 'U', 'G': 'G', 'C': 'C', 
                # Modified U
                'PSU': 'U', '5MU': 'U', 'H2U': 'U', '4SU': 'U',
                # Modified G
                'G7M': 'G', 'QUO': 'G', '2MG': 'G', 'OMG': 'G',
                'YYG': 'G', '7MG': 'G', 'M2G': 'G', 'YG':  'G',
                # Modified C
                'OMC': 'C', '5MC': 'C',
                # Modified A
                '1MA': 'A',
                }
chain_str = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

def output_write(output, outtype, base, x, y, z, serial=0, chain=0, resseq=0):
    if outtype == 'xyz':
        output.write("{:s} {:10.3f} {:10.3f} {:10.3f}\n".format(base, x, y, z))
    else:
        c = chain_str[(chain-1) % len(chain_str)]
        output.write('ATOM  {serial:5d} {atomname:4s} {resname:3s} {chain:1s}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}\n'.format(
                     serial = serial, atomname = '  '+base+' ', resname = ' '+base+' ',
                     chain = c, resseq = resseq, x = x, y = y, z = z))

x = 0.0
y = 0.0
z = 0.0
n = 0
nt_pre = None
ins_pre = None
flg_first = True
flg_next_nt = False
serial = 0
resseq = 0
chain = 0

for line in open(sys.argv[1], 'r'):
    
    if line.startswith('ATOM') or line.startswith('HETATM'):

        Atom = line[12:16].strip()
        nt = int(line[22:26])
        ins = line[26:27]
        nucleotide = line[17:21]
        #print('Nucleotide =', nucleotide)
        base = nucleotide.strip()
        #print('base = ', base)

        if Atom[-1] == "'" and not Atom[0] =='H':
            x += float(line[30:38])
            y += float(line[38:46])
            z += float(line[46:54])
            n += 1

        if flg_first:
            chain += 1
            resseq = 0
            nt_pre = nt
            ins_pre = ins
            prevbase = base 
            flg_first = False
            continue

        if nt != nt_pre or ins != ins_pre:
            flg_next_nt = True 

    elif line.startswith('TER'):
        #print('Finished Reading Chain')
        flg_next_nt = True
        flg_first = True

    else:
        continue

    # The line contains another nucleotide or TER line.
    if flg_next_nt:

        if n != 0:   # This condition is needed to exclude HETATM for HOH etc.

            if n != 9:
                print('Warning: the number of sugar atoms is not 9.', nt, ins, n)
        
            avgx = x/float(n)
            avgy = y/float(n)
            avgz = z/float(n)

            serial += 1
            resseq += 1
        
            try:
                output_write(output, outtype, unified_base[prevbase], avgx, avgy, avgz, serial=serial, chain=chain, resseq=resseq)

            except KeyError:
                print("error - Base = ", prevbase) 
                sys.exit(2)

            if flg_first and outtype == 'pdb':
                output.write('TER\n')

        x = 0.0
        y = 0.0
        z = 0.0
        n = 0

        nt_pre = nt
        ins_pre = ins
        prevbase = base 
        flg_next_nt = False

# Only in case there is no TER line at the end
if n != 0:

    if n != 9:
        print('Warning: the number of sugar atoms is not 9.', nt, ins, n)

    avgx = x/float(n)
    avgy = y/float(n)
    avgz = z/float(n)

    serial += 1
    resseq += 1

    try:
        output_write(output, outtype, unified_base[prevbase], avgx, avgy, avgz, serial=serial, resseq=resseq)

    except KeyError:
        print("error - Base = ", prevbase) 
        sys.exit(2)

if outtype == 'pdb':
    output.write('END\n')

output.close()

# If xyz, add the number of particles on top.
if outtype == 'xyz':
    lines = []
    for l in open(sys.argv[2]):
        lines.append(l)

    output = open(sys.argv[2], 'w')
    output.write('%i\n\n' % len(lines))
    for l in lines:
        output.write(l)
    output.close()
