#!/usr/bin/env python3
import sys,os,textwrap

#SEQDIR='/projects/dcking@colostate.edu/support_data/ce11'
SEQDIR='/projects/dcking@colostate.edu/support_data/ce11/Nmasked'

SEQDIR_CONTENTS = os.listdir(SEQDIR)
seqs = {}
infile = sys.argv[1]
lflank = 0
rflank = 0

if len(sys.argv) > 2:
    lflank = int(sys.argv[2])
    rflank = lflank
    if len(sys.argv) > 3:
        rflank = int(sys.argv[3])

for line in open(infile):
    fields = line.strip().split()
    chrom, start, end = fields[0], int(fields[1]) - lflank, int(fields[2]) + rflank

    # OPEN SEQUENCE FILE IF NOT ALREADY OPEN
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
            fh.close()
            

    # EXTRACT AND PRINT SEQUENCE
    # prevent indexing from tripping on out of range values
    start = max(start, 0)
    end = min(end, len(seqs[chrom]['seq']))
    flankText = ''
    if lflank > 0 or rflank > 0:
        flankText = " (-%d,+%d)" % (lflank,rflank)
    print(">%s:%d-%d%s" % (chrom,start,end,flankText))
    LINEWIDTH = 150
    print("\n".join(textwrap.wrap( seqs[chrom]['seq'][start:end], LINEWIDTH)))
    

