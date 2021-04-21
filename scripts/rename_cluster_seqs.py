#!/usr/bin/env python3

import sys
from Bio import SeqIO

records = SeqIO.parse(sys.argv[1], 'fasta')

uc = [l.strip() for l in open(sys.argv[2], 'r') if l.startswith('C')]


def find_cluster(s, uc):
    for l in uc:
        if s in l:
            return l.split('\t')[1]

for r in records:
    items = r.id.split('|')
    rid = "{}|{}".format(items[0], items[1])
    cluster = find_cluster(r.id, uc)
    print('>Cluster{} {}\n{}\n'.format(cluster, rid, r.seq))
