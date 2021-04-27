#! /usr/bin/env bash

python2 quast.py \
	--ref-bam $2  \
	--eukaryote \
	--large \
	--rna-finding \
	--fragmented \
	--conserved-genes-finding \
	-o genome_report \
	-r $1 \
	--threads 16 \
	$1

