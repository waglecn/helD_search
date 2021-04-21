#!/usr/bin/env python3

import sys
from Bio import SeqIO
import os

genome = sys.argv[1]

in_aa = f'hits/{genome}.hits'
in_up = f'fa/{genome}.upstream'

hits = SeqIO.to_dict(SeqIO.parse(in_aa, 'fasta'))
raes = SeqIO.to_dict(SeqIO.parse(in_up, 'fasta'))

for k in hits.keys():
	i = k.split('|')[1]
	print(raes[i].format('fasta'))

