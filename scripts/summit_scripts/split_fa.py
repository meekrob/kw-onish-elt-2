#!/usr/bin/env python3
import sys

try:
    mfa_file = sys.argv[1]
    mfa = open(sys.argv[1])
    out_base = '.'.join(mfa_file.split('.')[:-1])
    out_ext = mfa_file.split('.')[-1]
    split_by = open(sys.argv[2])

except FileNotFoundError:
    print("A file in", sys.argv[1:], "was not found", file=sys.stderr)
    sys.exit(1)
except IndexError: 
    print("USAGE: %s sequences.mfasta split.txt" % sys.argv[0], file=sys.stderr)
    print("arg1: sequences.mfasta: a multi-fasta file", file=sys.stderr)
    print("arg2: split.txt: levels to split the sequences into.", file=sys.stderr)
    print("  The number of sequences in sequences.mfasta must be the same as the number of lines in split.txt", file=sys.stderr)
    print("  The first WS-delimited column in the file split.txt will be taken as the factor to group the split by. It will be the suffix 'fac' to individual output files named 'sequences_fac.mfasta'", file=sys.stderr)
    sys.exit(1)

except:
    print("Unexpected error:", sys.exc_info()[0], file=sys.stderr)
    sys.exit(1)

fas = []
current_fa = None
for line in mfa:
    if not line.strip():
        continue

    if line.startswith('>'):
        if current_fa is not None: 
            fas.append( current_fa ) 
        current_fa = { 'header': line.lstrip('>'), 'seq': ''}
        continue

    current_fa['seq'] += line.strip()
    
    
if current_fa is not None: 
    fas.append( current_fa ) 

print("Read", len(fas), "sequences")

outfas = {}
for fa, factor_line in zip(fas, split_by.readlines()):
    factor = factor_line.strip().split()[0]
    if factor not in outfas:
        outname = "%s_%s.%s" % (out_base,factor,out_ext)
        print("opening", outname, "for", factor, sys.stderr)
        outfas[factor] = open(outname, "a")
    print(">" + fa['header'], file=outfas[factor], end='')
    print(fa['seq'], file=outfas[factor])

for k,v in outfas.items():
    print("closing file for %s" % k, file=sys.stderr)
    v.close()
    
    
