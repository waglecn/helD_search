#!/usr/bin/env python3

import sys

hitfile = sys.argv[1]
raefile = sys.argv[2]

hits = [r.strip().split('\t') for r in open(hitfile, 'r').readlines()]


default = ['N/A' for i in range(16)]

def parse_raes(raefile):
    raes = [
        r.strip().split(' ') for r in open(raefile, 'r').readlines()
        if not r.startswith('#')
    ]
    output = []
    for r in raes:
        temp = [i for i in r if i is not '']
        assert(len(temp) == 16)
        # shitty parsing
        # target_name, accession
        # query_name, accession
        temp = temp[:4] + [
            int(temp[4]),  # hmm_from
            int(temp[5]),  # hmm_to
            int(temp[6]),  # ali_from
            int(temp[7]),  # ali_to
            int(temp[8]),  # env_from
            int(temp[9]),  # env_to
            int(temp[10]),  # sq_len
            temp[11],      # strand (+ / -)
            float(temp[12]),  # E-value
            float(temp[13]),  # score
            float(temp[14]),  # bias
            temp[15],       # description of target
        ]
        if temp[13] > 13.0:
            output.append(temp)
    return output


raes = parse_raes(raefile)

import os

head, tail = os.path.split(hitfile)
genome = tail[:-7]

for h in hits:
    print(genome, h[0])
    print(raes)

