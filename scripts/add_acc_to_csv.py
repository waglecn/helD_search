#!/usr/bin/env python3

import sys

nodes = open(sys.argv[1], 'r').read()
nodes = [n.split(',') for n in nodes.split('\n')[:-1]]
length = len(nodes[0])

accs = open(sys.argv[2], 'r').read()
accs = [a.split(',') for a in accs.replace('"', '').split('\n')]

clusters = open(sys.argv[3], 'r').read()
clusters = [
    c.split('\t') for c in clusters.split('\n') if not c.startswith('#')
]

clust = {}
for c in clusters[:-1]:
    items = c[0].split(' ')
    cid = int(
        [i for i in items if '>Cluster' in i][0].replace('>Cluster', '')
    )
    if cid not in clust:
        clust[cid] = []
    clust[cid].append(c)
for c in clust:
    print(c)
    for i in clust[c]:
        print("\t{}".format(i[-3]))
exit()



for n in nodes[1:]:
    if len(n) != length:
        print(n)
    # print(n[:5], n[5:])
    result = [a for a in accs if a[0] in n[2]]
    assert len(result) > 0
    new = n[:5] + result[0] + n[5:]
    print(",".join(new))
