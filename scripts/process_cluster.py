#!/usr/bin/env python3
import sys
import os

def process(infile):
    centroids = {}
    clusters = {}

    lines = [l.strip().split('\t') for l in open(infile, 'r')]
    for l in lines:
        c = int(l[1])
        if l[0] == 'S':
            if c not in centroids:
                centroids[c] = l[8]
            if c not in clusters:
                clusters[c] = []
            if l[8] not in clusters[c]:
                clusters[c].append(l[8])
        if l[0] == 'H':
            clusters[c].append(l[8])
    return centroids, clusters


def main():
    centroids, clusters = process(sys.argv[1])

    outfile = sys.argv[1][:-3] + '.txt'
    print(outfile, file=sys.stderr)

    with open(outfile, 'w') as outh:
        for c in clusters:
            outh.write(
                'Cluster {0}, Length {1}:\t{2}'.format(
                    c, len(clusters[c]), centroids[c]
                )
            )
            for i in clusters[c][1:]:
                outh.write("\t{}\n".format(i))
            outh.write('\n')
    extras = sys.argv[2:]
    if len(extras) > 0 and extras[0] == 'clusters':
        # do something
        from Bio import SeqIO
        seqs = SeqIO.to_dict(SeqIO.parse('clustered_hits/all.sorted.fasta', 'fasta'))
        for i in extras[1:]:
            print(f'cluster {i}', file=sys.stderr)
            assert int(i) in clusters.keys()
            for j in clusters[int(i)]:
                print(seqs[j.split(' ')[0]].format('fasta'))
    elif len(extras) > 0 and extras[0] == 'find':
        rae_clusters = {}
        for c in clusters:
            cluster_count = 0
            for s in clusters[c]:
                if not s.endswith('|NA'):
                    if str(c) not in rae_clusters:
                        rae_clusters[str(c)] = 0
                    rae_clusters[str(c)] += 1
        for r in rae_clusters:
            print('{}\t{}\t{}\t{:.2f}'.format(
                r,
                len(clusters[int(r)]),
                rae_clusters[str(r)],
                rae_clusters[str(r)] / len(clusters[int(r)])
            ))

if __name__ == '__main__':
    main()
