#!/usr/bin/env python3

import sys

from ete3 import NCBITaxa

ncbi = NCBITaxa()

def get_desired_ranks(taxid):
    desired_ranks = ['phylum', 'genus', 'species']
    try:
        lineage = ncbi.get_lineage(taxid)
    except ValueError:
        lineage = []
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict(
        (rank, taxid) for (taxid, rank) in lineage2ranks.items()
    )
    return {
        '{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for
        rank in desired_ranks
    }

def main(taxid):
    rank_ids = get_desired_ranks(taxid)
    names = []
    for i in rank_ids.values():
        if isinstance(i, int):
            names.append(ncbi.translate_to_names([i])[0])
        else:
            names.append(i)
    return names

if __name__ == '__main__':
    names = main(sys.argv[1])
    print(','.join(names))