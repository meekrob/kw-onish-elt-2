#!/usr/bin/env python3
import sys

filename = sys.argv[1]
key_field = sys.argv[2]

if len(sys.argv) > 3:
    outfile = open(sys.argv[3], 'w')
else:
    outfile = sys.stdout

infile = open(filename)

header = infile.readline().strip()
header_fields = header.split("\t")
NFields = len(header_fields)
print("Header fields are:", "\t\n".join(header_fields), file=sys.stderr)

index_of_key = header_fields.index( key_field )
print("index of key field (%s): %d" % (key_field, index_of_key), file=sys.stderr)

last_fields = None
in_join = False
for line in infile:
    this_fields = line.strip().split("\t")
    this_key = this_fields[ index_of_key ]
    if last_fields:
        last_key = last_fields[ index_of_key ]
        if this_key == last_key:
            for i in range(NFields):
                last_field = last_fields[i]
                this_field = this_fields[i]
                if last_field != this_field or i > 17: # i > 17 because there was one with two of the same blockSize
                    if last_field == 'NA':
                        last_fields[i] = this_field
                    elif this_field == 'NA':
                        pass
                    else:
                        last_fields[i] += ',' + this_field
            
        else:
            print("\t".join(last_fields), file=outfile)
            last_fields = this_fields
    else:
        last_fields = this_fields

print("\t".join(last_fields), file=outfile)
                        
    
    

