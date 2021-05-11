#!/usr/bin/env python
"""
This script renames the query and subect fileds in the blast output table to an
integers representing their position in the all.fasta file, so they can be run
through the col_table.py script.

This script will also optionally filter for E-value at the specified limit

"""


import sys
from Bio import SeqIO

records = SeqIO.parse('../clustered_hits/all.fasta', 'fasta')
records = [r.id for r in records]

ELIMIT = 1E-75

indata = [l.strip().split('\t') for l in open('all.out', 'r')]
for i in indata:
    temp = [
        str(records.index(i[0])),
        str(records.index(i[1]))
    ] + i[2:]
    if float(temp[-2]) < ELIMIT:
        print("\t".join(temp))
