#!/usr/bin/env python3
"""
This script processes sequences from the genome flat files

Depending on the second argument, either 'faa' or 'upstream', will
extract protein sequences or the upstream 500 bp nucleotide sequences
for each genome.
"""

import sys
from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial


def main():
    input_file = sys.argv[1]

    encoding = guess_type(input_file)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

    seqs, upstreams = load(input_file, _open)

    if sys.argv[2] == 'faa':
        for s in seqs:
            print(s)
    elif sys.argv[2] == 'upstream':
        for u in upstreams:
            print(u)
    else:
        pass


def load(input_file, _open):
    seqs = []
    upstreams = []
    with _open(input_file) as f:
        records = [r for r in SeqIO.parse(f, 'gb')]
        for i, r in enumerate(records):
            CDS = [f for f in r.features if f.type == 'CDS']
            for j, c in enumerate(CDS):
                seq, upstream = process(c, r)
                seqs.append('>{}-{}\n{}\n'.format(
                    i, j, seq
                ))
                if len(upstream) > 0:
                    upstreams.append('>{}-{}\n{}\n'.format(i, j, upstream))

    return seqs, upstreams


def process(c, r):
    seq = c.extract(r).seq.translate(table=11)
    upstream = None
    if c.strand == 1:
        upstream = r[c.location.start - 500: c.location.start].seq
        # print(c.location.start, c.location.end, c.location, c.qualifiers)
    elif c.strand == -1:
        upstream = r[c.location.end:c.location.end + 500].seq

    return seq, upstream


if __name__ == '__main__':
    main()
