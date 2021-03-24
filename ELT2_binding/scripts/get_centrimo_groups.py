#!/usr/bin/env python
import sys, re, math
import numpy as np
from html.parser import HTMLParser
from io import StringIO
import xml.etree.ElementTree as ET

def main():

    # read meme data to calculate the ICs
    motifs = {}
    for motif in parse_meme_file( 'combined.meme' ):
        motifs[ motif['IDs'][1] ] = motif
    #print(*motifs.keys())

    # read tomtom results 
    tomtom = parse_tomtom_file('meme_tomtom_out/tomtom.xml') # a dict keyed like; ('RYKGGCGS', 'MEME-7') 
    tomtom.update( parse_tomtom_file('dreme_tomtom_out/tomtom.xml') )
    
    centrimo_data = []
    centrimo_indexes = []
    header = ''

    # centrimo output table
    for line in open('centrimo_out/centrimo.tsv'):
        if line.startswith('db_index'): 
            header = line.strip().split()
            continue
        if line.startswith('#'): continue
        if line.strip() == '': continue

        fields = line.strip().split()
        if fields[0] == '': del fields[0]

        db_index, motif_id = int(fields[0]), fields[1]
        motif_meme = fields[1] + '-' + fields[2]
        centrimo_indexes.append( (db_index, motif_id ) )
        centrimo_data.append( fields )
        
    # load a meme-chip.html that has been opened in a browser and OuterHTML saved with developer tools (option-Command-I)
    parser = memeChipHtmlParser()
    parser.feed( open( 'meme-chip-loaded.html' ).read())

    header.insert(4, 'information_content')
    header.insert(0, 'centrimo_group')
    header.insert(3, 'meaningful_id')
    print(*header, sep="\t")

    last_group_i = None

    for group_i, centrimo_group in parser.centrimo_groups:
        #print(centrimo_group)
        for identifier in centrimo_group:
            db,motif_id = identifier.split('+')
            db = int(db.replace('db','') ) + 1
            #print( centrimo_indexes.index( (db, motif_id) ), end=' ' )
            which_one = centrimo_indexes.index( (db, motif_id) )
            centrimo_fields = centrimo_data[ which_one ]
            meme_id = centrimo_fields[1] + '-' + centrimo_fields[2]
            ic = 0
            if meme_id not in motifs:
                print(line, file=sys.stderr)
                print(motif_meme, "not found", file=sys.stderr)
                #sys.exit(1)
                ic = -1 
            else:
                ic = motifs[meme_id]['ic']

            centrimo_fields.insert(4, "%.3f" % ic)
            if group_i != last_group_i:
                centrimo_fields.insert(0, 'top_%d' % group_i)
            else:
                centrimo_fields.insert(0, group_i)

            meaningful_id = centrimo_fields[3]
            try:
                tomtom_hit = tomtom[ (centrimo_fields[2],centrimo_fields[3]) ] # ('TCTTATCAGT', 'MEME-1')
                meaningful_id = tomtom_hit['match_alt']
            except:
                pass

            centrimo_fields.insert(3, meaningful_id)

            print(*centrimo_fields, sep="\t")
            last_group_i = group_i

class memeChipHtmlParser(HTMLParser):
    def __init__(self):
        self.stack = []

        self.in_a = False
        self.link_attr = None
        self.in_motifbox = False
        self.in_col_known_motifs = False
        self.group_i = 0
        self.centrimo_groups = []

        # for each div.motifbox:
        self.centrimo_id = None # the CentriMo link under "Discovery/Enrichment Program". Assumed to always be there.
        self.centrimo_group_data = None # the link, if it exists at the bottom of the motifbox with text "CentriMo Group"

        HTMLParser.__init__(self)

    def attr_to_dict( attrs ):
        attr_dict = {}
        for k,v in attrs: attr_dict[k] = v
        return attr_dict

    def handle_starttag(self, tag, attrs):
        attr_dict = memeChipHtmlParser.attr_to_dict( attrs )
        try: 
            ['img', 'br', 'meta', 'hr'].index(tag)
        except ValueError:
            self.stack.append( (tag,attr_dict) )

        if tag == 'a': 
            self.in_a = True
            self.link_attr = attr_dict

        elif tag == 'div' and 'class' in attr_dict and attr_dict['class'] == 'motifbox':
            self.in_motifbox = True

        elif tag == 'td' and 'class' in attr_dict and attr_dict['class'] == 'col_known_motifs':
            self.in_col_known_motifs = True

    def end_motifbox(self):
        
        if self.centrimo_group_data is not None:
            self.centrimo_groups.append((self.group_i, self.centrimo_group_data))
            self.group_i += 1
        elif self.centrimo_id is not None:
            self.centrimo_groups.append((self.group_i, [self.centrimo_id]))
            self.group_i += 1
        else:
            #print("Error: motifbox has no centrimo or group id")
            pass # the first one triggers because it's the "Expand All Clusters/Collapse All Clusters" header

        self.centrimo_group_data = None
        self.centrimo_id = None
        self.in_motifbox = False
    

    def handle_endtag(self, tag):
        if self.stack[-1][0] == tag:
            attr_dict = self.stack[-1][1]
            if tag == 'div' and 'class' in attr_dict and attr_dict['class'] == 'motifbox':
                self.end_motifbox()

            elif tag == 'td' and 'class' in attr_dict and attr_dict['class'] == 'col_known_motifs':
                self.in_col_known_motifs = False

            self.stack.pop()
        else:
            print("%s is on top but reached endtag: %s" % (self.stack[-1][0], tag), file=sys.stderr) 
            sys.exit(1)

        if tag == 'a': 
            self.in_a = False
            self.link_attr = None

    def process_CentriMo_Group_link(href):
        url_args = href.split('?')[-1]
        group_args = url_args.split('show=')
        group_ids = []
        for group_arg in group_args:
            if not group_arg: continue
            group_id = group_arg.replace('&','')
            group_ids.append( group_id )

        return group_ids

    def handle_data(self, data):
        if self.in_motifbox:
            if self.in_a: 
                if self.in_col_known_motifs:
                    pass
                elif data.startswith("CentriMo Group"):
                    self.centrimo_group_data = memeChipHtmlParser.process_CentriMo_Group_link(self.link_attr['href'])
                        
                elif data.strip() == 'CentriMo':                            # <a title="Click to view the CentriMo output" 
                    self.centrimo_id = self.link_attr['href'].split('show=')[1]  # href="centrimo_out/centrimo.html?show=db0+SGCGCGMCGC">
                                                                            # CentriMo</a>
                

def parse_tomtom_file(xmlfile, e_threshold = .05):

    result_dict = {}

    tree = ET.parse(xmlfile)
    root = tree.getroot()
    format=root.tag

    row_keys = ['id', 'alt', 'tool_evalue','match_id', 'match_alt', 
                'tomtom_eval', 'tomtom_qval'] 

    queries = []
    for m in root.findall('queries/motif'):
        #print(m.get('id'), m.get('alt'))
        queries.append( m )
         #<motif db="0" id="AARYGCGC" alt="DREME-6" length="8" nsites="173" evalue="2.3e-018">
    
    targets = []
    for m in root.findall('targets/motif'):
        targets.append( m )
        #print(m.get('id'), m.get('alt'))

    for q in root.findall('matches/query'):
        query_idx = int(q.get('idx'))
        query = queries[query_idx]
        evalue = np.float( query.get('evalue') )
        if evalue > e_threshold:
            continue

        for hit in list(q):
            hit_idx = int(hit.get('idx'))
            target = targets[ hit_idx ]
            evalue = float(hit.get('ev'))
            qvalue = float(hit.get('qv'))

            row = [query.get('id'), query.get('alt'), query.get('evalue'),target.get('id'), target.get('alt'), hit.get('ev'), hit.get('qv')]
            rowd = {}
            for k,v in zip(row_keys, row):
                rowd[k] = v
            overall_key = (rowd['id'], rowd['alt'])
            result_dict[overall_key] = rowd
            break

    return result_dict
    
# end def parse_tomtom_file()
                    

#c_elegans_bg = { 'A': .25, 'C': .25, 'G': .25, 'T': .25 }
c_elegans_bg = { 'A': 3.246e-01, 'C': 1.754e-01, 'G': 1.754e-01, 'T': 3.246e-01 }

                         # symbol for revcomp
twofold = { 'AC' : 'M',  # TG: K
            'AG' : 'R',  # CT: Y
            'AT' : 'W',  # TA: W
            'CG' : 'S',  # CG: S
            'CT' : 'Y',  # GA: R
            'GT' : 'K' } # CA: M

revcomp_table = str.maketrans('ACGTMRWSYK', 'TGCAKYWSRM')


def print_motif(motif):
    print(*['MOTIF'] + motif['IDs'])
    print()

    print('letter-probability matrix:', end=' ')
    for k,v in motif['lpm_data'].items():
        print('%s= %s' % (k,v), end=' ')
    print('ic= %.3f' % float(motif['ic'])) # meme2meme will delete this field 
    for row in motif['lpm']:
        for b in [ 'A', 'C', 'G', 'T' ]:
            if b == 'T': 
                end='\n'
            else:
                end=' '
            print(" %.6f" % row[b], end=end)
    print()

    if 'URL' in motif: print(motif['URL'])
    print()
        

def read_lpm(fhi):
    mx_txt = ''
    while True:
        line = next(fhi)
        if not line.strip(): break
        mx_txt += line

    c = StringIO(mx_txt)
    lpm = np.loadtxt(c, dtype={'names': ('A', 'C', 'G', 'T'), 'formats': ('f4', 'f4', 'f4', 'f4')})
    return lpm

def get_ic(fg, bg):
    if fg == 0:
        if bg != .25: # must apply pseudocounts, assuring fg > 0
            raise ValueError
        return 0
    else:
        return fg * math.log( fg/bg, 2 )
    
def information_content(lpm):
    ic = 0
    for i in range( len(lpm) ):
        row_ic = 0
        for b in [ 'A', 'C', 'G', 'T' ]:
            row_ic += get_ic( lpm[b][i], c_elegans_bg[b] )
            
        ic += row_ic
    return ic

def parse_letter_probability_matrix_line(lpm_label):
    txt_fields = re.split(' |\=', lpm_label) # spaces and equal signs
    fields = filter(lambda x: x, txt_fields) # delete empty fields
    attr = {}
    for name in fields: # they now alternate btw name and value
        value = next(fields)
        attr[name] = value

    return attr

def apply_pseudocount(lpm, nsites, pseudocount=1):
    new_nsites = nsites+(4*pseudocount)
    for row in lpm:
        for b in [ 'A', 'C', 'G', 'T' ]:
            row[b] = ((row[b] * nsites)+pseudocount) / new_nsites

    return lpm, new_nsites
    

def parse_meme_file(filename):
    # the print statements (may be commented out) in here print data that precedes the motifs
    in_meme = filename
    fhi = open( in_meme )

    motif = None
    motifs = []
    ic_fails = 0
    while True:
        try:
            line = next(fhi)
        except StopIteration:
            break

        if line.startswith('MOTIF'): 
            if motif: motifs.append( motif )
            motif = {'IDs' : line.strip().split()[1:] }

        elif line.startswith('letter-probability matrix'): 
            lpm_label = line.strip().split(':')[1]

            motif['lpm_data'] = parse_letter_probability_matrix_line( lpm_label )

            lpm = read_lpm(fhi)
            try:
                ic = information_content(lpm)
            except ValueError: # raised if any cell in lpm is 0 and bg is not uniform (.25)
                #print(*["0 encountered in LPM but bg is not uniform. Applying pseudocounts to"] + motif['IDs'], file=sys.stderr)
                # apply pseudocount
                nsites = int(motif['lpm_data']['nsites'])
                lpm, new_nsites = apply_pseudocount(lpm, nsites, nsites/100)
                motif['lpm_data']['nsites'] = str(new_nsites)
                ic = information_content(lpm)
            motif['lpm'] = lpm
            motif['ic'] = ic

        elif line.startswith('Background letter frequencies'):
            #print(line, end='')
            #print(next(fhi))
            next(fhi)
        elif line.startswith('URL'):
            motif['URL'] = line.strip()
        elif line.strip():
            #print(line)
            pass

        if fhi.closed: break
            
    if motif: motifs.append( motif )

    return motifs

def revcomp(seq):
    return seq.translate(revcomp_table)[::-1]

def get_two_fold(pair):
    s = ''.join( sorted( pair ) )
    return twofold[s]

def consensus_char(row_tuples):
    row_tuples.sort(key=lambda r: -r[0])
    
    if row_tuples[0][0] >= .5: 
        return row_tuples[0][1]
    elif row_tuples[0][0] + row_tuples[1][0] >= .75:
        return get_two_fold(row_tuples[0][1] + row_tuples[1][1])
    else:
        return 'N'

def consense(pwm):
    lines = [ "A\tC\tG\tT" ]
    cons = ''
    for i in range(len(pwm)):
        row_tuples = list(zip(pwm[i], ['A','C','G','T'])) 
        c = consensus_char( row_tuples )
        lines.append( "%.4f\t%.4f\t%.4f\t%.4f\t%s" % (pwm[i][0], pwm[i][1], pwm[i][2], pwm[i][3], c))
        cons += c
    return cons,lines

if __name__ == '__main__':
    main()
