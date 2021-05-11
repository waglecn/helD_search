#!/usr/bin/env python3
"""
This script generates the fraction of rae statistic for each cluster, and adds
the value to the end of each tsv row

input: tsv file, ie clustered_hits/50/detected_helD.tsv

output: stdout tsv text
"""

import sys
import subprocess

infile = sys.argv[1]


def fraction(c):
    fasta_file = "detected_helD.c{}.fasta".format(c)
    r = [l.strip() for l in open(fasta_file, 'r')]
    total = len([l for l in r if l.startswith('>')])
    non_rae = len([l for l in r if l.endswith('|NA')])

    return (1 - (float(non_rae) / float(total)))


records = [r.strip().split('\t') for r in open(infile, 'r')]
for r in records:
    c = r[0]
    r += ["{:.3f}".format(fraction(c))]
    print('\t'.join(r))
