rule all:
    input:


rule trim_short:
    input:
        F = "reads/short/{short}.F.fq.gz",
        R = "reads/short/{short}.R.fq.gz"
    output:
        F = "reads/short_trimmed/{short}.R1.fq.gz",
        R = "reads/short_trimmed/{short}.R2.fq.gz"
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
        length_req = "--length_required 50",   #TODO FIGURE THIS VALUE OUT
        q = "-q 15",
        u = "-u 50",
        correction = "--correction",
    shell:
        """
        fastp --thread {threads} --in1 {input.F} --in2 {input.R} -h {log.html} -j {log.json} {params} &> {log.txt}
        """


rule kmergenie_short:
    input:
        F = "reads/short/{short}.R1.fq.gz",
        R = "reads/short/{short}.R2.fq.gz"
    output:
    #TODO figure out output
        "kmergenie/something"
    params:

    message:
        """
        Finding optimal Kmers for short reads with Kmergenie
        """
    shell:
        """
        kmergenie {input}
        """


rule kmergenie_long:
    input:
        "reads/long/{long}.fq.gz"
    output:

    params:

    message:
        """
        Finding optimal Kmers for long reads with Kmergenie
        """
    shell:
        """
        kmergenie {input}
        """


rule sparseassembler:
    input:
        in1 = "reads/short_trimmed/{short}.R1.fq.gz",
        in2 = "reads/short_trimmed/{short}.R2.fq.gz",
        kmer = "{kmergenie_short.output}"
    output:
        "sparseassembler/contigs.txt"
    params:
        in1 = "i1 ../{input.in1}",
        in2 = "i2 ../{input.in2}",
        genome_size = "GS 1780000000",
        node_thresh = "NodeCovTh 2",
        edge_thresh = "EdgeCovTh 1",
        g_thresh = "g 15",
        ld_val = "LD 0",
    message:
        """
        Assembling short reads with Kmergenie-derived K
        """
    shell:
        """
        mkdir sparseassembler && cd sparseassembler
        KMER=$(grep ../{kmer} somevalue)
        SparseAssembler k $KMER {params}
        """


rule dbg2olc:
    input:
        sparse = "sparseassembler/contigs.txt",
        longreads = "reads/long/{long}.fq.gz",
        kmergenie = "{kmergenie_long.output}"
    output:
        "dbg2olc/backbone_raw.fasta"
    message:
        """
        Assembling long and short reads with DBG2OLC
        """
    params:
        ld_val = "LD 0",
        adaptive_theta = "AdaptiveTh 0.001",
        kmer_cov_thresh = "KmerCovTh 3",
        min_overlap = "MinOverlap 50",
        rm_chimera = "RemoveChimera 1",
        infile = "Contigs ../{input.sparse}",
        longreads = "f ../{input.longreads}"
    shell:
        """
        mkdir dbg2olc && cd dbg2olc
        KMER=$(grep ../{kmergenie} somevalue)
        DBG2OLC k $KVAL {params}
        """
