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

RAE

GGGGCTTGCGGCAAGGCCC


- systematic search of RAE vs resistance +- HelR
- find other palindromes
- directional RAE HMM
venezuela locus tag 5092 - palindrome upstream
CCCCGACTACAATTGACAGTCGGGG

venezuelae
GCF_008639165.1_ASM863916v1_genomic