#! /usr/bin/env bash

python software/quast/quast.py \
	--ref-bam genome_report/YFT_genome.bam  \
	--eukaryote \
	--large \
	--rna-finding \
	--fragmented \
	--conserved-genes-finding \
	-o genome_report \
	-r polish/YFT.genome.fasta \
	--threads 16 \
	polish/YFT.genome.fasta

