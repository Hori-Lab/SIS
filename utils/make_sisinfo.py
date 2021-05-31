#!/usr/bin/env python

nchains = 64
nrepeat = 47
nmp_chain = 3 * nrepeat
nmp = nmp_chain * nchains

iserial = 0

fout = open('CAG47.sisinfo', 'w')

ibd = 0
for ichain in range(nchains):

    for imp in range(1, nmp_chain):

        ibd += 1
        imp_total = nmp_chain * ichain + imp
        fout.write('bond {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:f} {:f}\n'.format(ibd, imp_total, imp_total+1, ichain+1, imp, ichain+1, imp+1, 15.0, 5.9))

iangl = 0
for ichain in range(nchains):

    for imp in range(1, nmp_chain-1):

        iangl += 1
        imp_total = nmp_chain * ichain + imp
        fout.write('angl {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:d} {:f} {:f}\n'.format(iangl, imp_total, imp_total+1, imp_total+2, ichain+1, imp, ichain+1, imp+1, ichain+1, imp+2, 10.0, 2.618))

fout.close()

