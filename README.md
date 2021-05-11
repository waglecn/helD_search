# helD_search

This project is a set of scripts that will be used to search for homologues
of helD, a helicase proteien from B. subtilis being investigated for links to
resistance to rifamycin.

see https://www.uniprot.org/uniprot/O32215 for the original sequences
https://www.uniprot.org/uniprot/A0QUE0
https://www.biorxiv.org/content/10.1101/2020.07.20.211821v1.full.pdf

Requirements:
- usearch (in path) (www.drive5.com)


```
git clone https://github.com/waglecn/helD_search.git
```

```
conda env create -n helD -f environment.yaml
conda activate helD
```

RAE sequence from Spannogiannopolous et al, see resources rae.txt

GGGGCTTGCGGCAAGGCCC

Usage:

snakemake -c use-conda COMMAND

where COMMANDS are one or more of:
  
  download - downloads the gbff genome files into a common location. This path will need to be modified in the snakefile and other scripts.

  prep - uses usearch to cluster the sequences obtained from the B. subtilis and M. smegmatis HelD BLAST searches in /resources. Prepares a redundant set of sequences to search against the RefSeq genomes.

    This step also creates the .faa and .upstream files for each genome in the RefSeq list

  search - performs the identification of the RAE sequences in the .upstream data for each genome. Also performs the similarity search using the 99% redundant M. smegmatis and B. subtilis query sequence data against the RefSeq proteins.
  
  all_hits - performs all the downstream processing on the RefSeq helD sequences. This includes clustering, and alignment and tree building. Depending on the output, these are the clusters which may be used for CDD searchs and are the base files from which input tables to Cytoscape may be produced.

  summary - produces the csv summary file for each genome.


TODO:
- find other palindromes (not complete)
- directional RAE HMM (not complete)
S. venezuelae locus tag 5092 - palindrome upstream
CCCCGACTACAATTGACAGTCGGGG

Notes:
Example HMMER output for genome with 3 known RAE sequences
GCF_008639165.1_ASM863916v1_genomic
Query:       rae  [M=19]
Scores for complete hits:
    E-value  score  bias  Sequence  start    end  Description
    ------- ------ -----  --------  -----  -----  -----------
    0.00065   21.6   3.6  0-461       130    148  
    0.00065   21.6   3.6  0-462       423    441  
       0.01   17.9   3.5  0-3382      423    439  
  ------ inclusion threshold ------
        1.4   11.2   2.3  0-1425      311    323  
        2.2   10.6   1.1  0-906       454    465  
        3.2   10.1   3.2  0-5936       79     96  
          4    9.8   1.9  0-6442       16     27  
        4.4    9.7   2.0  0-6440      356    367 