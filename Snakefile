rule all:
    input: "misassembly/YFT_to_purgeI.bam"


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
        cut_qual = "--cut_mean_quality 15",
        adapters = "--detect_adapter_for_pe",
        length_req = "--length_required 100",
        q = "-q 15",
        u = "-u 50",
        correction = "--correction"
    shell:
        """
        fastp --thread {threads} --in1 {input.F} --in2 {input.R} --out1 {output.F} --out2 {output.R} -h {log.html} -j {log.json} {params} &> {log.txt}
        """


rule kmergenie_short:
    input:
        F = "reads/short_trimmed/{short}.illumina.R1.fq",
        R = "reads/short_trimmed/{short}.illumina.R2.fq"
    output:
        best_k = "kmergenie/short/{short}.best.k"
    log:
        k_report = "kmergenie/short/{short}_report.html"
    params:
        kmin = "-l 20",
        kmax = "-k 120",
        out_prefix = "-o {short}"
    threads: 16
    message:
        """
        Finding optimal Kmers for short reads with Kmergenie
        """
    shell:
        """
        mkdir -p kmergenie/short/ && cd kmergenie/short
        echo -e "../../{input.F}\n../../{input.R}" > shortreads.txt
        ../../software/kmergenie-1.7051/kmergenie shortreads.txt -t {threads} {params} > ../../{output.best_k}
        """


rule kmergenie_long:
    input:
        long_reads = "reads/long/{long}.pb.fasta",
        #short_contigs = "sparseassembler/{long}_contigs.txt"
    output:
        best_k = "kmergenie/long/{long}.best.k"
    log:
        k_report = "kmergenie/long/{long}_report.html"
    params:
        kmin = "-l 20",
        kmax = "-k 200",
        out_prefix = "-o {long}"
    threads: 16  
    message:
        """
        Finding optimal Kmers for long reads + sparseassembler contigs with Kmergenie
        """
    shell:
        """
        mkdir -p kmergenie/long/ && cd kmergenie/long/
        echo -e "../../{input.long_reads}" > ../../kmergenie/long/longreads.txt
        ../../software/kmergenie-1.7051/kmergenie longreads.txt -t {threads} {params} > ../../{output.best_k}
        """

rule sparseassembler:
    input:
        in1 = "reads/short_trimmed/{short}.illumina.R1.fq",
        in2 = "reads/short_trimmed/{short}.illumina.R2.fq",
        kmer = "kmergenie/short/{short}.best.k"
    output:
        contigs = "sparseassembler/{short}_contigs.txt"
    params:
        genome_size = "GS 2000000000",
        #node_thresh = "NodeCovTh 2",
        edge_thresh = "EdgeCovTh 1",
        g_thresh = "g 15",
        ld_val = "LD 0"
    message:
        """
        Assembling short reads with Kmergenie-derived K
        """
    shell:
        """
        KMER=$(grep "^best k:" {input.kmer} | grep -o '[^ ]*$')
        KCOV=$(grep "for best k:" {input.kmer} | grep -o '[^ ]*$')
        mkdir -p sparseassembler
        cd sparseassembler
        SparseAssembler k $KMER NodeCovTh $KCOV i1 ../{input.in1} i2 ../{input.in2} {params}
        mv Contigs.txt ../{output}
        """


rule dbg2olc:
    input:
        sparse = "sparseassembler/{prefix}_contigs.txt",
        longreads = "reads/long/{prefix}.pb.fasta",
        kmer = "kmergenie/long/{prefix}.best.k"
    output:
        contigs = "dbg2olc/{prefix}_backbone_raw.fasta",
        contig_info = "dbg2olc/{prefix}_consensus_info.txt"
    message:
        """
        Assembling long and short reads with DBG2OLC
        """
    params:
        ld_val = "LD 0",
        kmer = "k 31",
        k_coverage = "KmerCovTh 2",
        adaptive_theta = "AdaptiveTh 0.01",
        min_overlap = "MinOverlap 50",
        rm_chimera = "RemoveChimera 1"
    shell:
        """
        mkdir -p dbg2olc
        cd dbg2olc
        #KMER=$(grep "^best k:" ../{input.kmer} | grep -o '[^ ]*$')
        #KCOV=$(grep "for best k:" ../{input.kmer} | grep -o '[^ ]*$')
        DBG2OLC Contigs ../{input.sparse} f ../{input.longreads} {params} > dbg.log
        mv backbone_raw.fasta ../{output.contigs}
        mv DBG2OLC_Consensus_info.txt ../{output.contig_info}
        """

rule concat_contigs:
    input:
        longreads = "reads/long/{prefix}.pb.fasta",
        short_contigs = "sparseassembler/{prefix}_contigs.txt",
    output:
        concat_contigs = "reads/consensus_concat/{prefix}_contigs.pb.fasta"
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
        Using BLASR + Sparc to perform a consensus
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
        ../software/dbg2olc/split_and_run_sparc.sh $DBG_CONT $CONT_INF $CONTIGS $TMPDIR 2 {threads} > ../{log}
        mv tmp/final_assembly.fasta ../{output.consensus} && rm -r tmp
        """

rule map_for_purge:
    input:
        in1 = "reads/short_trimmed/{prefix}.illumina.R1.fq",
        in2 = "reads/short_trimmed/{prefix}.illumina.R2.fq",
        consensus = "consensus/{prefix}_consensus.fasta"
    output: 
        mapfile = "purge_haplotigs/first/{prefix}_to_consensus.bam",
        mapindex = "purge_haplotigs/first/{prefix}_to_consensus.bam.bai"
    params: 
        samfile = "purge_haplotigs/first/{prefix}_to_consensus.sam"
    message: "Mapping short reads onto the consensus genome"
    threads: 16
    shell:
        """
        software/bwa-mem2/bwa-mem2 index {input.consensus}
        software/bwa-mem2/bwa-mem2 mem -t {threads} {input.consensus} {input.in1} {input.in2} > {params.samfile}
        software/bwa-mem2/sam2bam {params.samfile} {threads}
        """

rule purge_haplotigs_I_hist:
    input:
        consensus = "consensus/{prefix}_consensus.fasta",
        mapfile = "purge_haplotigs/first/{prefix}_to_consensus.bam"
    output:
        histo = "purge_haplotigs/first/{prefix}_to_consensus.bam.gencov",
        hist_image = "purge_haplotigs/first/{prefix}_to_consensus.bam.histogram.png"
    message: "Generating coverage histogram for purging"
    threads: 16
    params:
        depth = "-d 620"
    shell:
        """
        cd purge_haplotigs/first/
        purge_haplotigs hist -b ../../{input.mapfile} -g ../../{input.consensus} -t {threads} {params.depth}
        """

rule purge_haplotigs_suspects_I:
    input:
        hist_cov = "purge_haplotigs/first/{prefix}_to_consensus.bam.gencov"  
    output:
        cov_out = "purge_haplotigs/first/{prefix}_coverage_stats.csv"
    params:
        low = "-low 192",
        mid = "-mid 352",
        high = "-high 448"
    message: "Finding suspect contigs"
    shell:
        """
        purge_haplotigs cov -i {input} {params} -o {output} 
        """

rule purge_haplotigs_I:
    input:
        consensus = "consensus/{prefix}_consensus.fasta",
        mapfile = "purge_haplotigs/first/{prefix}_to_consensus.bam",
        suspects = "purge_haplotigs/first/{prefix}_coverage_stats.csv"
    output:
        curated = "purge_haplotigs/first/{prefix}_purge_I.fasta"
    log:
        haplotigs = "purge_haplotigs/first/{prefix}_purge_I.haplotigs.fasta",
        artefacts = "purge_haplotigs/first/{prefix}_purge_I.artefacts.fasta",
        reassignments = "purge_haplotigs/first/{prefix}_purge_I.reassignments.tsv",
        logs = "purge_haplotigs/first/{prefix}_purge_I.contig_associations.log"
    threads: 16
    message: "Purging haplotigs"
    params:
        prefix = "-o {prefix}_purge_I"
    shell:
        """
        cd purge_haplotigs/first
        purge_haplotigs purge {params} -t {threads} -g ../../{input.consensus} -c ../../{input.suspects} -d -b ../../{input.mapfile}
        """

rule map_for_MEC:
    input:
        in1 = "reads/short_trimmed/{prefix}.illumina.R1.fq",
        in2 = "reads/short_trimmed/{prefix}.illumina.R2.fq",
        purged_asm = "purge_haplotigs/first/{prefix}_purge_I.fasta"
    output: 
        mapfile = "misassembly/{prefix}_to_purgeI.bam",
        mapindex = "misassembly/{prefix}_to_purgeI.bam.bai"
    params: 
        samfile = "misassembly/{prefix}_to_purgeI.sam"
    message: "Mapping short reads onto the purged genome for misassembly removal"
    threads: 16
    shell:
        """
        software/bwa-mem2/bwa-mem2 index {input.purged_asm}
        software/bwa-mem2/bwa-mem2 mem -t {threads} {input.purged_asm} {input.in1} {input.in2}
        #software/bwa-mem2/sam2bam {params.samfile} {threads}
        #software/BamQC/bin/bamqc {output.mapfile} -t {threads} -q -o misassembly
        """