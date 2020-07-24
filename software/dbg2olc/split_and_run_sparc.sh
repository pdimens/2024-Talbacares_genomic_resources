#!/usr/bin/env bash

###USAGE
#split_and_run_sparc.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR] [ITERATIONS] [THREADS]

backbone_fasta=$1
consensus_fasta=$2
reads_fasta=$3
split_dir=$4
iterations=$5
THREADS=$6

mkdir -p $split_dir

SPLITREADS=$(realpath ../software/dbg2olc/split_reads_by_backbone_openclose.py)
python2 $SPLITREADS -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} 

for file in $(find ${split_dir} -name "*.reads.fasta"); do
    chunk=`basename $file .reads.fasta`

    for iter in `seq 1 ${iterations}`; do
        blasr --nproc 16 ${split_dir}/${chunk}.reads.fasta ${split_dir}/${chunk}.fasta --bestn 1 -m 5 --minMatch 19 --out ${split_dir}/${chunk}.mapped.m5 &> /dev/null
        Sparc m ${split_dir}/${chunk}.mapped.m5 b ${split_dir}/${chunk}.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o ${split_dir}/${chunk} 
    done
    
    #remove after use to save space
    rm ${split_dir}/${chunk}.mapped.m5 ${split_dir}/${chunk}.reads.fasta
    cat ${split_dir}/${chunk}.consensus.fasta >> ${split_dir}/final_assembly.fasta && rm ${split_dir}/${chunk}.consensus.fasta
done