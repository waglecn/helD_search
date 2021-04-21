#!/usr/bin/env python3
import sys
import ete3
from Bio import SeqIO

names_file = sys.argv[1]
tree_file = sys.argv[2]

seqs = SeqIO.to_dict(SeqIO.parse(names_file, 'fasta'))
names = list(seqs.keys())


tree = ete3.Tree(tree_file)

for l in tree.get_leaves():
    for k in names:
        if k.startswith(l.name):
            l.name = seqs[k].description
            # print(l.name, seqs[k].description)


tree.write(outfile='detected_helD.trimal.fasttree.rename.tree')