#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 2014/06/17
@author: Naoto Hori
'''

import sys
import argparse
from dcd import DcdFile
from sisout import SisoutFile
from sisbp import SisbpFile
#import resource

NUM_OUTFILE_OPEN = 100

#resource.setrlimit(resource.RLIMIT_NOFILE, (2000,2000))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Convert REMD simulation result from files for each replica to files for each label.',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('indir', help='Input directory path')
    parser.add_argument('name', help='Filename prefix, e.g. "md" for md_0001.out, md_0002.out ...')
    parser.add_argument('nrep', type=int, help='Number of replicas')
    parser.add_argument('outdir',help='Output directory path')
    parser.add_argument('--bp', action="store_true", help='Flag for bp file. Use this when you need to convert .bp file as well.')

    args = parser.parse_args()

    #Open input files
    in_sisout_files = []
    in_dcd_files = []
    in_sisbp_files = []
    for id_rep in range(1, args.nrep+1):
        sisout = SisoutFile('%s/%s_%04i.out' % (args.indir, args.name, id_rep))
        sisout.open_to_read()
        #sisout.read_header()   ## Commented out since the header will be read later
        in_sisout_files.append(sisout)

        dcd = DcdFile('%s/%s_%04i.dcd' % (args.indir, args.name, id_rep))
        dcd.open_to_read()
        #dcd.read_header()  ## Commented out since the header will be read later
        in_dcd_files.append(dcd)

        if args.bp:
            bp = SisbpFile('%s/%s_%04i.bp' % (args.indir, args.name, id_rep))
            bp.open_to_read()
            #dcd.read_header()  ## Commented out since the header will be read later
            in_sisbp_files.append(bp)

    id_finish = 0
    while id_finish < args.nrep:

        id_begin_now = id_finish + 1
        id_end_now = id_finish + NUM_OUTFILE_OPEN
        if id_end_now > args.nrep:
            id_end_now = args.nrep
        id_finish = id_end_now

        # Prepare output files
        list_id_out = list(range(id_begin_now, id_end_now+1))
        out_sisout_files = {}
        out_dcd_files = {}
        out_sisbp_files = {}
        for id_lab in list_id_out:
            out = SisoutFile('%s/%s_%04i.out' % (args.outdir, args.name, id_lab))
            out.open_to_write()
            out_sisout_files[id_lab] = out
            dcdout = DcdFile('%s/%s_%04i.dcd' % (args.outdir, args.name, id_lab))
            dcdout.open_to_write()
            out_dcd_files[id_lab] = dcdout
            if args.bp:
                bpout = SisbpFile('%s/%s_%04i.bp' % (args.outdir, args.name, id_lab))
                bpout.open_to_write()
                out_sisbp_files[id_lab] = bpout


        # Rewind (by reading header again)
        for id_rep in range(1, args.nrep+1):
            idx = id_rep - 1

            in_sisout_files[idx].read_header()
            in_dcd_files[idx].read_header()
            if args.bp:
                in_sisbp_files[idx].read_header()
 
            if (id_rep in list_id_out):
                out_sisout_files[id_rep].copy_header(in_sisout_files[idx])
                out_sisout_files[id_rep].write_header()
                out_dcd_files[id_rep].set_header( in_dcd_files[idx].get_header() )
                out_dcd_files[id_rep].write_header()
                if args.bp:
                    out_sisbp_files[id_rep].copy_header( in_sisbp_files[idx] )
                    out_sisbp_files[id_rep].write_header()

        while in_dcd_files[0].has_more_data():
            for idx in range(args.nrep):
                (outdata, outlines) = in_sisout_files[idx].read_onestep()
                coord_matrix = in_dcd_files[idx].read_onestep()
                if args.bp:
                    pairs, energies = in_sisbp_files[idx].read_onestep()

                label = int( outdata[0][in_sisout_files[idx].head_col.label] )

                if label in list_id_out:
                    out_sisout_files[label].write_onestep( outlines )
                    out_dcd_files[label].write_onestep( coord_matrix )
                    if args.bp:
                        out_sisbp_files[label].write_onestep( pairs, energies )

        # Close output files
        for f in list(out_sisout_files.values()):
            f.close()
        for f in list(out_dcd_files.values()):
            f.close()
        if args.bp:
            for f in list(out_sisbp_files.values()):
                f.close()


    # Close input files
    for f in in_sisout_files:
        f.close()
    for f in in_dcd_files:
        f.close()
    if args.bp:
        for f in in_sisbp_files:
            f.close()
