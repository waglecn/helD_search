import os
import glob

"""
-Goal: create a phylogenetic tree of the RIF resistance AAA Helicases related
to a sequence found in M. smegmatis.

Problem: These helicases are large, and multiple homologues exist in Gram+ve
organisms.

Solution: Search for ill-defined homlogues. HelD belongs to super family 1
(SF1) - which mainly has homlogy in the middle part of the protein. This may
be due to ATPase, so it will be difficult to identify all homologues.

Using HelD (M.smeg) as query, get top 5000 hits with evalue 1e-255 (0.0,
effectively), then cluster these into families. The centroid of each family
will be used to search a set of RefSeq genomes from Firmicutes and
Actinobacteria for the closest homologue over a series of similarity
thresholds (ct, cluste threshold).

One particular ct will best cluster these families, but how to tell? look for
association with RAE, as derived from Spannogiannopolous et al converted into
a nucleotide HMM.

This HMM will be queried against the upstream region of all putative HelD hits.
The cluster threshold which combines the greateset number of RAE+ HelD
seqeunces will then be used to identify the sequences to be used for the MSA
and then the phylogenetic tree.

"""


# configfile: "config.yaml"
taxa = ['Actinobacteria', 'Firmicutes']
org = ['msme', 'bsub']

cluster_thresholds = [99, 95, 93, 90, 85, 80, 75, 70, 65, 60, 55, 50]
genomes = glob.glob("/media/nick/2TB/backup-ignore/ngd/refseq/bacteria/*.gbff.gz")
# genomes = genomes[:10]
genomes = [os.path.split(g[:-8])[-1] for g in genomes]


rule all:

rule prep:
    input:
        # Cluster Bsub and Msme hits
        expand(
            'usearch/{org}-{input_seqs}.{ct}.centroids.fasta',
            input_seqs='HelD-5000', ct=cluster_thresholds,
            org=['msme', 'bsub']
        ),
        # prep target genomes, amino acid sequences and upstream regions
        expand(
            'fa/{genome}.faa', genome=genomes
        ),
        expand('fa/{genome}.upstream', genome=genomes),

rule search:
    # identify Bsub and Msme helD homologues and rae sequences
    input:
        expand("fa/{genome}.nhmmer.txt", genome=genomes),
        expand(
            'search/{org}/99/{genome}.search.txt', genome=genomes, org=org
        ),

rule summary:
    # summarize the HelD hits
    input:
        expand("csv/{genome}.csv", genome=genomes),

rule all_hits:
    # process the hits
    input:
        "clustered_hits/all.sorted.fasta",
        expand(
            "clustered_hits/all.{ct}.sorted.centroids", ct=cluster_thresholds
        ),
        expand(
            "clustered_hits/all.{ct}.sorted.txt", ct=cluster_thresholds
        ),
        expand("clustered_hits/{ct}/detected_helD.tsv", ct=[50, 55, 60, 65, 70]),
        expand("clustered_hits/{ct}/detected_helD.all.trimal.fasttree.tree", ct=[50, 55, 60, 65, 70])


rule download:
    input:
        expand('download/{taxon}.target.taxlist.meta', taxon=taxa),


###############################################################################
# External Data Rules
###############################################################################
rule download_refseq_taxlist:
    threads: 1
    output:
        "download/{taxon}.target.taxlist"
    shell:
        "python scripts/gimme-taxa.py -j "
        " {wildcards.taxon} > {output}"

rule update_download_store:
    threads: 3
    input:
        "download/{taxon}.target.taxlist"
    params:
        out_location = "/media/nick/2TB/backup-ignore/ngd/refseq/"
    output:
        "download/{taxon}.target.taxlist.meta"
    log:
        "download/{taxon}.target.taxlist.meta.log"
    shell:
        "ngd -v -s refseq -r3 -p 3 -o {params.out_location} "
        "-t {input} --flat-output -F genbank "
        "-m {output} -l complete bacteria 2>&1 | tee {log}"


##############################################################################
# genome rules
##############################################################################

rule process_fa_genome:
    threads: 1
    input:
        "/media/nick/2TB/backup-ignore/ngd/refseq/bacteria/{genome}.gbff.gz"
    output:
        "fa/{genome}.faa"
    shell:
        "scripts/extract_fa.py {input} faa > {output}"

rule process_upstream_genome:
    threads: 1
    input:
        "/media/nick/2TB/backup-ignore/ngd/refseq/bacteria/{genome}.gbff.gz"
    output:
        "fa/{genome}.upstream"
    shell:
        "scripts/extract_fa.py {input} upstream > {output}"


##############################################################################
# input cluster rule
##############################################################################

rule length_sort:
    threads: 1
    input:
        "resources/{seqs}.fasta"
    output:
        "resources/{seqs}.sorted.fasta"
    shell:
        "scripts/usearch -sortbylength {input} -fastaout {output}"


rule cluster_input:
    threads: 1
    input:
        "resources/{seqs}.sorted.fasta"
    output:
        centroids = "usearch/{seqs}.{ct}.centroids.fasta",
        uc = "usearch/{seqs}.{ct}.uc"
    shell:
        "scripts/usearch -cluster_fast {input} -id 0.{wildcards.ct} "
        "-centroids {output.centroids} -uc {output.uc} "

##############################################################################
# BLAST rules
##############################################################################

rule make_blast_db:
    threads: 1
    input:
        "usearch/{org}-HelD-5000.{ct}.centroids.fasta"
    output:
        "usearch/{org}-HelD-5000.{ct}.centroids.fasta.dmnd"
    shell:
        "diamond makedb -d {input} --in {input}"

rule search_db:
    threads: 1
    input:
        db = "usearch/{org}-HelD-5000.{ct}.centroids.fasta.dmnd",
        query = "fa/{genome}.faa"
    output:
        "search/{org}/{ct}/{genome}.search.txt"
    shell:
        # "blastp -db usearch/HelD-5000.{wildcards.ct}.centroids.fasta "
        # "-query {input.query} -outfmt 6 -evalue 1e-255 -max_target_seqs 1 "
        # "-num_threads {threads} > {output}"
        "diamond blastp -p {threads} --max-target-seqs 1 "
        "-e 1e-100 -d {input.db} "
        "-q {input.query} --outfmt 6  > {output}"

##############################################################################
# nhmmer rules
##############################################################################

rule search_rae:
    threads: 1
    input:
        "fa/{genome}.upstream"
    output:
        "fa/{genome}.nhmmer.txt"
    log:
        "fa/{genome}.nhmmer.log"
    shell:
        "nhmmer --cpu {threads} --toponly --incT 13.0 --tblout {output} "
        "resources/rae.fasta {input} 2>&1 | tee {log} "

##############################################################################
# summarize by genome
##############################################################################
rule summarize_hit:
    threads: 1
    output:
        "csv/{genome}.csv"
    shell:
        'mkdir -p hits ;'
        "scripts/per_genome.py {wildcards.genome} > {output}"

rule combine_hits:
    threads: 1
    output:
        "clustered_hits/all.fasta"
    shell:
        "cat hits/*.hits > {output}"

rule hits_length_sort:
    threads: 1
    input:
        "clustered_hits/all.fasta"
    output:
        "clustered_hits/all.sorted.fasta"
    shell:
        "scripts/usearch -sortbylength {input} -fastaout {output}"

rule cluster_helD_hits:
    threads: 1
    input:
        "clustered_hits/all.sorted.fasta"
    output:
        centroids = "clustered_hits/all.{ct}.sorted.centroids",
        uc = "clustered_hits/all.{ct}.sorted.uc"
    shell:
        "scripts/usearch -cluster_fast {input} -id 0.{wildcards.ct} "
        "-centroids {output.centroids} -uc {output.uc} -target_cov 0.75"

rule process_uc:
    threads: 1
    input:
        "clustered_hits/all.{ct}.sorted.uc"
    output:
        "clustered_hits/all.{ct}.sorted.txt"
    shell:
        "scripts/process_cluster.py {input}"

rule summarize_clusters:
    threads: 1
    input:
        "clustered_hits/all.{ct}.sorted.uc"
    output:
        all = "clustered_hits/{ct}/detected_helD.all.fasta",
        tsv = "clustered_hits/{ct}/detected_helD.tsv"
    shell:
        "scripts/process_cluster.py {input} find > {output.tsv} ; "
        "for i in $(cut -f 1  {output.tsv}) ; do "
            "scripts/process_cluster.py {input} clusters $i > "
            "clustered_hits/{wildcards.ct}/detected_helD.c$i.fasta ; "
        " done ;"
        "cat clustered_hits/{wildcards.ct}/detected_helD.c*.fasta > {output.all}"

rule align_clusters:
    threads: 8
    input:
        uc = "clustered_hits/all.{ct}.sorted.uc",
        tsv = "clustered_hits/{ct}/detected_helD.50.tsv"
    shell:
        "for i in $(cut -f 1  {intput.tsv}) ; do "
            "mafft --auto --thread {threads} "
            "clustered_hits/{wildcards.ct}/detected_helD.c$i.fasta > "
            "clustered_hits/{wildcards.ct}/detected_helD.c$i.mafft.align ; "
        " done ;"


rule align_ct_all:
    threads: 8
    input:
        "clustered_hits/{ct}/detected_helD.all.fasta"
    output:
        "clustered_hits/{ct}/detected_helD.all.mafft.align"
    shell:
        "mafft --auto --thread {threads} {input} > {output}"

rule trimal_all:
    threads: 8
    input:
        "clustered_hits/{ct}/detected_helD.all.mafft.align"
    output:
        "clustered_hits/{ct}/detected_helD.all.trimal.align"
    shell:
        "trimal -automated1 -in {input} -out {output}"

rule fasttree_all:
    threads: 8
    input:
        "clustered_hits/{ct}/detected_helD.all.trimal.align"
    output:
        "clustered_hits/{ct}/detected_helD.all.trimal.fasttree.tree"
    shell:
        "fasttreeMP -wag < {input} > {output}"

