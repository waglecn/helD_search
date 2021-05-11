#!/usr/bin/env python3

import sys
from Bio import SeqIO
from progress.bar import Bar
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


infasta = sys.argv[1]
records = [r for r in SeqIO.parse(infasta, 'fasta')]
records_trunc = [r.id for r in records]
records_full = [r.description for r in records]

print(records_trunc[0])
print(records_full[0])

inblast = sys.argv[2]
blast = [l.strip().split('\t') for l in open(inblast, 'r')]
blast = [b for b in blast if b[0] != b[1]]
blast = [b for b in blast if float(b[-1]) > 700.0]

scores = [float(b[-1]) for b in blast]
# with open('scores.out', 'w') as outh:
#     for s in scores:
#         outh.write("{s}\n".format(s))
sns.histplot(scores)
plt.show()

pairs = set()
pbar = Bar('Processing blast', max=len(blast))
for b in blast:
    p_0 = records_trunc.index(b[0])
    p_1 = records_trunc.index(b[1])
    p = min((p_1, p_0), (p_0, p_1))
    pairs.add(p)
    pbar.next()
pbar.finish()
print(len(blast), len(pairs))

incluster = sys.argv[3]
uc = [l.strip().split('\t') for l in open(incluster, 'r')]
uc_centroids = [records_full.index(u[-2]) for u in uc if u[0] == 'C']
cluster_map = {}
for u in  uc:
    if u[0] in ['H', 'S']:
        cluster_map[u[-2]] = u[1]

inraeclust = sys.argv[4]
raeclust = [r.strip().split("\t") for r in open(inraeclust, 'r')]

nodes = {}
for p in pairs:
    if p[0] not in nodes:
        nodes[p[0]] = []
    if p[1] not in nodes:
        nodes[p[1]] = []

# fill in nodes table
columns = [
    'id', 'full_name', 'assembly', 'coords', 'phylum', 'taxon', 'rae',
    'cluster', 'clust_type'
]

ctypes = {
    '8': 'Gtf-HelD-fusion',
    '12': 'Gtf-HelD-fusion',
    '53': 'HelD',
    '68': 'HelD',
    '74': 'HelD',
    '124': 'HelD',
    '133': 'HelD',
    '195': 'HelD',
    '227': 'HelD_partial',
    '230': 'HelD',
    '235': 'HelD',
    '239': 'HelD_partial',
    '243': 'HelD_partial',
    '247': 'HelD_partial',
    '258': 'HelD_partial',
    '259': 'HelD_partial',
    '266': 'HelD',
    '270': 'HelD_partial',
    '280': 'HelD_partial',
    '313': 'HelD_partial',
    '319': 'HelD_partial',
    '342': 'HelD_partial',
    '346': 'Gtf',
    '354': 'HelD_partial',
    '361': 'Gtf',
    '367': 'Gtf',
    '370': 'Gtf',
    '381': 'HelD_partial',
    '384': 'HelD_partial',
    '391': 'HelD_partial',
    '405': 'HelD_partial',
    '414': 'HelD_partial',
    '417': 'HelD_partial',
}
node_table = []
for n in nodes.keys():
    full_name = records_full[n]

    split_name = full_name.split('|')
    cluster = 'NA'
    if full_name in cluster_map:
        cluster = cluster_map[full_name]
    ctype = 'None'
    if cluster in ctypes:
        ctype = ctypes[cluster]

    node_table.append(
        [n, full_name] + split_name + [cluster, ctype]
    )

with open('pairs.out', 'w') as outh:
    for p in pairs:
        outh.write("{}\t{}\n".format(p[0], p[1]))

table = pd.DataFrame(node_table)
table.columns = columns
table.to_csv('nodes.csv')
