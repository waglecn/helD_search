#!/usr/bin/env python3

import sys
from Bio import SeqIO
from ete3 import NCBITaxa

ncbi = NCBITaxa()

infile = sys.argv[1]
print(infile)
records = [r for r in SeqIO.parse(infile, 'fasta')]
names = [r.description.split('|')[-2] for r in records]
name_dict = ncbi.get_name_translator(names)
for n in names:
    while True:
        try:
            print(name_dict[n][0])
            break
        except KeyError:
            new = n.split(' ')[:-1]
            print(f'cannot find \'{n}\', propose \'{new}\'', file=sys.stderr)
            if len(new) > 0:
                n = ' '.join(new)
                new_names = ncbi.get_name_translator([n])
                name_dict = {**name_dict, **new_names}
                print(f'trying \'{n}\'', file=sys.stderr)
            else: 
                break


