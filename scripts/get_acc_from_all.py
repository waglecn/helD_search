#!/usr/bin/env python3
import sys
import pathlib
from Bio import SeqIO
import gzip

genomes = {}
records = SeqIO.parse(sys.argv[1], 'fasta')


def process_genome(genome):
    gfile = "{}.gbff.gz".format(genome)
    gpath = pathlib.Path(
        '/media/nick/2TB/backup-ignore/ngd/refseq/bacteria/', gfile
    )
    inh = gzip.open(gpath,'rt')
    records = SeqIO.parse(inh, 'genbank')
    contigs = []
    for r in records:
        contigs.append([f for f in r.features if f.type == 'CDS'])
    return contigs



for r in records:
    items = r.id.split('|')
    genome = items[0]
    ca = items[1].split('-')
    c = int(ca[0])
    a = int(ca[1])
    if genome not in genomes:
        genomes[genome] = []
    genomes[genome].append((c, a))


for g in genomes:
    contigs = process_genome(g)
    for item in genomes[g]:
        # print(contigs[item[0]][item[1]])
        locus_tag = contigs[item[0]][item[1]].qualifiers['locus_tag'][0]
        note = ""
        try:
            protein_id = contigs[item[0]][item[1]].qualifiers['protein_id'][0]
        except KeyError:
            protein_id = ""
            note = contigs[item[0]][item[1]].qualifiers['note'][0]

        print("\"{}\",\"{}\",\"{}\",\"{}\"".format(
            "{}|{}-{}".format(g, item[0], item[1]),
            locus_tag,
            protein_id,
            note
        ))


