rule all:
	input: "init_polish/YFT.racon_longshort.fasta"
	message: "Generating report for genome assembly"
	threads: 16
	#shell: "python software/quast/quast.py --ref-bam {input.mapfile} --eukaryote --large --rna-finding --conserved-genes-finding -o polish -r {input} --threads {threads} {input}"


rule trim_short:
	input:
		F = "reads/short/{short}.illumina.F.fq.gz",
		R = "reads/short/{short}.illumina.R.fq.gz"
	output:
		F = "reads/short_trimmed/{short}.illumina.R1.fq",
		R = "reads/short_trimmed/{short}.illumina.R2.fq"
	log:
		json = "reads/short_trimmed/{short}.json",
		html = "reads/short_trimmed/{short}.html",
		txt = "reads/short_trimmed/{short}.trim.log"
	message:
		"""
		Trimming raw short reads with fastp
		"""
	threads: 16
	params:
		cut_type = "--cut_front --cut_tail", 
		cut_window = "--cut_window_size 5",
		cut_qual = "--cut_mean_quality 20",
		adapters = "--detect_adapter_for_pe",
		length_req = "--length_required 50",
		q = "-q 15",
		u = "-u 50",
		correction = "--correction"
	shell:
		"""
		fastp --thread {threads} --in1 {input.F} --in2 {input.R} --out1 {output.F} --out2 {output.R} -h {log.html} -j {log.json} {params} &> {log.txt}
		"""

rule sparseassembler:
	input:
		in1 = "reads/short_trimmed/{short}.illumina.R1.fq",
		in2 = "reads/short_trimmed/{short}.illumina.R2.fq",
	output:
		contigs = "sparseassembler/{short}_contigs.txt"
	params:
		genome_size = "GS 2000000000",
		node_thresh = "NodeCovTh 2",
		edge_thresh = "EdgeCovTh 1",
		g_thresh = "g 15",
		k = "k 51",
		chims = "ChimeraTh 2",
		chims_thresh = "ContigTh 2"
	message:
		"""
		Assembling short reads with {params.k}
		"""
	shell:
		"""
		mkdir -p sparseassembler
		cd sparseassembler
		SparseAssembler i1 ../{input.in1} i2 ../{input.in2} {params}
		mv Contigs.txt ../{output}
		"""

rule dbg2olc:
	input:
		sparse = "sparseassembler/{prefix}_contigs.txt",
		longreads = "reads/long/{prefix}.pb.fasta",
	output:
		contigs = "dbg2olc/{prefix}_backbone_raw.fasta",
		contig_info = "dbg2olc/{prefix}_consensus_info.txt"
	message:
		"""
		Assembling long and short reads with DBG2OLC
		"""
	params:
		ld_val = "LD 0",
		kmer = "k 17",
		k_coverage = "KmerCovTh 2",
		adaptive_theta = "AdaptiveTh 0.001",
		min_overlap = "MinOverlap 15",
		rm_chimera = "RemoveChimera 1"
	shell:
		"""
		mkdir -p dbg2olc
		cd dbg2olc
		DBG2OLC Contigs ../{input.sparse} f ../{input.longreads} {params} > dbg.log
		mv backbone_raw.fasta ../{output.contigs}
		mv DBG2OLC_Consensus_info.txt ../{output.contig_info}
		"""

rule concat_contigs:
	input:
		longreads = "reads/long/{prefix}.pb.fasta",
		short_contigs = "sparseassembler/{prefix}_contigs.txt",
	output:
		concat_contigs = temp("reads/consensus_concat/{prefix}_contigs.pb.fasta")
	message: "Concatenating contigs for consensus"
	shell:
		"""
		cat {input} > {output}
		"""

rule consensus:
	input:
		dbg_contigs = "dbg2olc/{prefix}_backbone_raw.fasta",
		contig_info = "dbg2olc/{prefix}_consensus_info.txt",
		concat_contigs = "reads/consensus_concat/{prefix}_contigs.pb.fasta"
	output:
		consensus = "consensus/{prefix}_consensus.fasta",
	log: 
		"consensus/{prefix}_consensus.log"
	message:
		"""
		Using BLASR + pbdagcon to perform a consensus
		"""
	threads: 16
	shell:
		"""
		mkdir -p consensus/tmp
		DBG_CONT=$(realpath {input.dbg_contigs})
		CONT_INF=$(realpath {input.contig_info})
		CONTIGS=$(realpath  {input.concat_contigs})
		TMPDIR=$(realpath consensus/tmp)
		cd consensus
		../software/dbg2olc/split_and_run_pbdagcon.sh $DBG_CONT $CONT_INF $CONTIGS $TMPDIR 2 {threads} > ../{log}
		mv tmp/final_assembly.fasta ../{output.consensus} && rm -r tmp
		"""

rule map_long_racon_polish:
	input:
		reads = "reads/long/{prefix}.pb.fasta",
		consensus = "consensus/{prefix}_consensus.fasta"
	output: 
		mapfile = temp("init_polish/{prefix}_to_consensus.sam"),
	message: "Mapping long reads onto the consensus genome for first round of racon polishing"
	threads: 10
	shell:
		"""
		minimap2 -t {threads} -ax map-pb --sam-hit-only {input.consensus} {input.reads} > {output}
		"""

rule racon_polish_long:
	input:
		reads = "reads/long/{prefix}.pb.fasta",
		assembly = "consensus/{prefix}_consensus.fasta",
		mapfile = "init_polish/{prefix}_to_consensus.sam"
	output:
		asm = "init_polish/{prefix}.racon_long.fasta"
	message: "Polishing consensus with racon using long reads"
	threads: 10
	shell:
		"""
		racon -t {threads} {input.reads} {input.mapfile} {input.assembly} > {output}
		"""

rule subsample_shortreads:
	input:
		reads_f = "reads/short_trimmed/{prefix}.illumina.R1.fq",
		reads_r = "reads/short_trimmed/{prefix}.illumina.R2.fq",
	output:
		reads_f = "reads/short_trimmed/{prefix}.illumina.70x.R1.fq",
		reads_r = "reads/short_trimmed/{prefix}.illumina.70x.R2.fq",
	message: "Subsampling 25% of the short reads with seqtk"
	shell:
		"""
		seqtk sample -s100 {input.reads_f} .25 > {output.reads_f}
	seqtk sample -s100 {input.reads_f} .25 > {output.reads_f}
		"""

rule racon_preprocess_illumina:
	input:
		reads_f = "reads/short_trimmed/{prefix}.illumina.70x.R1.fq",
		reads_r = "reads/short_trimmed/{prefix}.illumina.70x.R2.fq",
	output: temp("reads/short_trimmed/{prefix}.illumina.70x.racon.fq")
	message: "Preprocessing short reads to work for racon polishing"
	shell: "software/racon/racon_preprocess.py {input} > {output}"

rule map_short_racon_polish:
	input:
		reads = "reads/short_trimmed/{prefix}.illumina.70x.racon.fq",
		assembly = "init_polish/{prefix}.racon_long.fasta"
	output:
		mapfile = temp("init_polish/{prefix}.short_to_longpolished.sam")
	message: "Mapping 70x subsampled short reads to long-read polished consensus"
	threads: 10
	shell:
		"""
		minimap2 -t {threads} -ax sr --sam-hit-only {input.assembly} {input.reads} | samtools view -h -F4 -q10 -@{threads} > {output}
		"""

rule racon_short_polish:
	input:
		reads = "reads/short_trimmed/{prefix}.illumina.70x.racon.fq",
		mapfile = "init_polish/{prefix}.short_to_longpolished.sam",
		assembly = "init_polish/{prefix}.racon_long.fasta"
	output:
		asm = "init_polish/{prefix}.racon_longshort.fasta"
	message: "Polishing consensus with racon using subsampled short reads"
	threads: 10
	shell:
		"""
		racon -t {threads} {input.reads} {input.mapfile} {input.assembly} > {output}
		"""
