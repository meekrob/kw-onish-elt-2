#!/usr/bin/env python3
import sys,os,textwrap

SEQDIR='/projects/dcking@colostate.edu/support_data/ce11'
SEQDIR_CONTENTS = os.listdir(SEQDIR)
seqs = {}
infile = sys.argv[1]
for line in open(infile):
    fields = line.strip().split()
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])

    if chrom not in seqs :
        if chrom + ".fa" not in SEQDIR_CONTENTS:
            print("warning: no seqfile for chromosome", chrom, file = sys.stderr)
        else:
            filename = os.path.join(SEQDIR, chrom + ".fa")
            print("opening %s to read seq" % filename, file=sys.stderr)
            fh = open(filename)
            header = fh.readline().strip()
            seq = fh.read().replace("\n", "")
            seqs[chrom] = { 'header': header, 'seq': seq }
            

    print("> %s:%d-%d" % (chrom,start,end))
    print("\n".join(textwrap.wrap( seqs[chrom]['seq'][start:end], 150)))
    

