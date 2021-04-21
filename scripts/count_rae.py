#!/usr/bin/env python3
import os
import glob
import statistics
from per_genome import parse_raes
from tax_csv import main as get_names
import subprocess

files = glob.glob("./fa/*.txt")

all_output = []
counts = []
non_counts = []
phyla = []

Act_frac_rae = []

for f in files:
    o = parse_raes(f)
    counts.append(len(o))
    if len(o) > 0:
        non_counts.append(len(o))
    for i in o:
        all_output.append(i)

    genome = os.path.split(f)[-1].replace('.nhmmer.txt', '')
    cmd = [
        './scripts/get_taxid.sh',
        '/media/nick/2TB/backup-ignore/ngd/refseq/bacteria/{}.gbff.gz'.format(
            genome
        )
    ]
    taxid = subprocess.check_output(cmd).decode('utf-8').strip()
    names = get_names(taxid)
    phyla.append(names)
    if names[0] == 'Actinobacteria':
        Act_frac_rae.append(len(o))

num_act = len([i for i in phyla if i[0] == 'Actinobacteria'])
num_fir = len([i for i in phyla if i[0] == 'Firmicutes'])
num_other = len(phyla) - num_act - num_fir


print(
    len(all_output),
    statistics.mean(counts),
    len(non_counts),
    statistics.mean(non_counts),
    num_act, num_fir, num_other,
    sum(Act_frac_rae)
)


