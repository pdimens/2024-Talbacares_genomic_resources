rule all:
    input: "dbg2olc/YFT_backbone_raw.fasta"


rule trim_short:
    input:
        F = "reads/short/{short}.illumina.F.fq.gz",
        R = "reads/short/{short}.illumina.R.fq.gz"
    output:
        F = "reads/short_trimmed/{short}.illumina.R1.fq.gz",
        R = "reads/short_trimmed/{short}.illumina.R2.fq.gz"
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
        length_req = "--length_required 50",
        q = "-q 15",
        u = "-u 50",
        correction = "--correction"
    shell:
        """
        fastp --thread {threads} --in1 {input.F} --in2 {input.R} --out1 {output.F} --out2 {output.R} -h {log.html} -j {log.json} {params} &> {log.txt}
        """


rule kmergenie_short:
    input:
        F = "reads/short/{short}.illumina.R1.fq.gz",
        R = "reads/short/{short}.illumina.R2.fq.gz"
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
        echo -e "{input.F}\n{input.R}" > shortreads.txt
        ../../software/kmergenie-1.7051/kmergenie shortreads.txt -t {threads} {params} > {output.best.k} 
        """


rule kmergenie_long:
    input:
        long_reads = "reads/long/{long}.pb.fasta",
        short_contigs = "sparseassembler/{long}_contigs.txt"
    output:
        best_k = "kmergenie/long/{long}.best.k"
    log:
        k_report = "kmergenie/long/{long}_report.html"
    params:
        kmin = "-l 30",
        kmax = "-k 200",
        out_prefix = "-o {long}"
    threads: 16  
    message:
        """
        Finding optimal Kmers for long reads + sparseassembler contigs with Kmergenie
        """
    shell:
        """
        mkdir -p kmergenie/long/ && cd kmergenie/long
        echo -e "{input.long_reads}\n{input.short_contigs}" > longreads.txt
        ../../software/kmergenie-1.7051/kmergenie longreads.txt -t {threads} {params} > {output.best.k} 
        """

rule decompress:
    input:
        in1 = "reads/short_trimmed/{short}.illumina.R1.fq.gz",
        in2 = "reads/short_trimmed/{short}.illumina.R2.fq.gz"
    output:
        in1 = "reads/short_trimmed/{short}.illumina.R1.fq.",
        in2 = "reads/short_trimmed/{short}.illumina.R2.fq"
    message:
        """
        Decompressing short reads for SparseAssembler
        """
    shell:
        """
        gunzip {input}
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
        node_thresh = "NodeCovTh 2",
        edge_thresh = "EdgeCovTh 1",
        g_thresh = "g 15",
        ld_val = "LD 0"
    message:
        """
        Assembling short reads with Kmergenie-derived K
        """
    shell:
        """
        mkdir sparseassembler && cd sparseassembler
        KMER=$(grep "^best k:" {kmer} | grep -o '[^ ]*$')
        SparseAssembler k $KMER i1 ../{input.in1} i2 ../{input.in2} {params}
        mv Contigs.txt {short}_contigs.txt
        """


rule dbg2olc:
    input:
        sparse = "sparseassembler/{prefix}_contigs.txt",
        longreads = "reads/long/{prefix}.pb.fasta",
        kmergenie = "kmergenie/long/{prefix}.best.k"
    output:
        contigs = "dbg2olc/{prefix}_backbone_raw.fasta"
    message:
        """
        Assembling long and short reads with DBG2OLC
        """
    params:
        ld_val = "LD 0",
        adaptive_theta = "AdaptiveTh 0.001",
        #kmer_cov_thresh = "KmerCovTh 3",
        min_overlap = "MinOverlap 50",
        rm_chimera = "RemoveChimera 1"
    shell:
        """
        mkdir dbg2olc && cd dbg2olc
        KMER=$(grep "^best k:" {kmer} | grep -o '[^ ]*$')
        KCOV=$(grep "for best k:" {kmer} | grep -o '[^ ]*$')
        DBG2OLC k $KMER KmerCovTh $KCOV Contigs ../{input.sparse} f ../{input.longreads} {params}
        mv backbone_raw.fasta {prefix}_backbone_raw.fasta
        """
