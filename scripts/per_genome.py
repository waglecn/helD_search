#!/usr/bin/env python3
"""
This script generates the sequence file with the updated sequence definition
lines annotated with taxonomy and RAE information

"""

import sys
import os
import subprocess
from Bio import SeqIO
from tax_csv import main as get_names


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
        if temp[13] > 10.0 and (temp[5] - temp[4]) >= 14.0:
            output.append(temp)
    return output


def parse_rae_align(raealignfile):
    aligns = [
        r.strip()for r in open(raealignfile, 'r').readlines()
        if not r.startswith('#')
    ]
    targets = [r[2:].strip() for r in aligns if r.startswith('>>')]
    temp = {}
    for t in targets:
        items = [l.strip().split(' ') for l in aligns if l.startswith(t)][0]
        temp[t] = [i for i in items if i is not '']

    return temp


def main():
    genome = sys.argv[1]
    m_blast_hits_file = "search/msme/99/{}.search.txt".format(genome)
    b_blast_hits_file = "search/bsub/99/{}.search.txt".format(genome)
    nhmmer_table = "fa/{}.nhmmer.txt".format(genome)
    nhmmer_out = 'fa/{}.nhmmer.log'.format(genome)

    for f in [
        m_blast_hits_file, b_blast_hits_file, nhmmer_table, nhmmer_out
    ]:
        assert os.path.exists(f), print(f)

    records = SeqIO.to_dict(SeqIO.parse('./fa/{}.faa'.format(genome), 'fasta'))

    b_hit = [l.strip().split('\t') for l in open(b_blast_hits_file, 'r')]
    m_hit = [l.strip().split('\t') for l in open(m_blast_hits_file, 'r')]

    hit = list(set([h[0] for h in b_hit + m_hit]))

    raes = parse_raes(nhmmer_table)

    aligns = parse_rae_align(nhmmer_out)

    cmd = [
        './scripts/get_taxid.sh',
        '/media/nick/2TB/backup-ignore/ngd/refseq/bacteria/{}.gbff.gz'.format(
            genome
        )
    ]
    taxid = subprocess.check_output(cmd).decode('utf-8').strip()
    names = get_names(taxid)

    row = {}
    fasta_out = ""
    for h in hit:
        query = ''
        if h in [h[0] for h in b_hit]:
            query += 'B'
        if h in [h[0] for h in m_hit]:
            query += 'M'
        # from, to, score, align
        rae = ['NA', 'NA', 'NA', 'NA']
        r = [r for r in raes if h == r[0]]
        if len(r) > 0:
            r = r[0]
            temp = [r[4], r[5], r[13], None]

            if h in aligns:
                temp[3] = aligns[h][2]
                # print(temp, aligns[h[0]])
                temp[3] = '{}{}{}'.format(
                    'n' * (1 - temp[0]),
                    temp[3],
                    'n' * (19 - temp[1])
                )
            rae = temp

        temp = [genome, taxid] + names + [query, h] + rae
        row[h] = temp
        print(','.join([str(t) for t in temp]))
        fasta_out += ">{}|{}|{}|{}|{}\n{}\n".format(
            genome,
            h,
            temp[2], temp[4], temp[-1],
            records[h].seq)
    if len(fasta_out) > 0:
        with open('hits/{}.hits'.format(genome), 'w') as outh:
            outh.write(
                fasta_out.replace('<not', 'not').replace('present>', 'present')
            )


if __name__ == '__main__':
    main()
