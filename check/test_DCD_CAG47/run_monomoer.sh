./monomer.py -f "(CAG)47" \
             -x 100 -n 10000 \
             -t md.dcd -o md.out \
             -H 1.67 \
             -i CAG47.pdb \
             -K 150.0 --cutoff 50.0 \
             1> out 2> err
