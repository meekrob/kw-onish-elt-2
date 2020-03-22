#!/usr/bin/env python3
import sys

colors = {
  'protein_coding': '0,0,0',
  'pseudogene': '128,128,128',
  'tRNA': '166,206,227',
  'miRNA': '31,120,180',
  'ncRNA': '178,223,138',
  'rRNA': '51,160,44',
  'snRNA': '251,154,153',
  'snoRNA': '227,26,28',
  'lincRNA': '253,191,111',
  'piRNA': '255,127,0',
  'antisense_RNA': '202,178,214'}

def safe_map_int(fields):
    out = []
    for f in fields:
        try:
            out.append(int(f))
        except ValueError:
            pass
    return f

for line in open( sys.argv[1] ):
    fields = line.strip().split("\t")

    fields[5] = min(map(int,fields[5].split(',')))
    fields[6] = max(map(int,fields[6].split(',')))
    fields[7] = len(fields[7].split(','))

    fields.insert(4, "1000")
    fields.insert(8, "0,0,0")
    # 19 and 20 go to 10 and 11
    fields.insert(10, fields[-2])
    fields.insert(11, fields[-1])

    # Need to check the block starts/ends against the given chrom start/end
    # wormbase allows different transcripts to be different sections within
    # the (probably) comprehensive range of the gene.
    # UCSC requires the chrom start/end to match the range covering the block 
    # starts/ends
    chromStart = int(fields[1])
    chromEnd = int(fields[2])
    thickStart = int(fields[6])
    thickEnd = int(fields[7])
    blockSizes = list(map(int, fields[10].split(',')))
    blockStarts = list(map(int, fields[11].split(',')))
    if fields[5] == '-':
        blockSizes.reverse()
        blockStarts.reverse()
    
    if blockStarts[0] != 0:
        chromStart = chromStart + blockStarts[0]
        blockStarts[0] = 0

    chromEnd = chromStart + blockStarts[-1] + blockSizes[-1] 
        
    fields[1] = str(chromStart)
    fields[2] = str(chromEnd)
    if thickStart == thickEnd:
        fields[6] = fields[1]
        fields[7] = fields[1]
    else:
        fields[6] = str( max(chromStart,thickStart) )
        fields[7] = str( min(chromEnd,thickEnd) )

    fields[10] = ",".join( map(str,blockSizes) )
    fields[11] = ",".join( map(str,blockStarts) )

    # use the wormbase locus ID as the name, I have put the wbps_transcript_id in place when not set
    fields[3] = fields[18]
    # color code the biotype
    fields[8] = colors[ fields[17] ]

    print("\t".join(map(str, fields[:-2])))
