#!/bin/bash
#module load fastqc

sample_prefix=$1
FASTQDIR="/projects/ps-epigen/seqdata/"
OUTPUT_DIR="/projects/ps-epigen/outputs/libQCs/${sample_prefix}"

mkdir -p $OUTPUT_DIR

for i in `seq 1 2`
do
    # do fastqc on raw fastq files

    pre="${sample_prefix}_R${i}"
    [[ ! -f "${FASTQDIR}${pre}.fastq.gz" ]] &&  pre=${pre}.trim
    pf="${FASTQDIR}${pre}.fastq.gz"
    echo "running fastqc $pf ..."        
    fastqc  -t 16 -o $OUTPUT_DIR $pf

    # do fastq_screen on trimmed fastq files (perfer)

    echo "running fastscreen $p ..."
    pre="${sample_prefix}_R${i}.trim"
    [[ ! -f "${FASTQDIR}${pre}.fastq.gz" ]] &&  pre=${pre/.trim/}
    pf="${FASTQDIR}${pre}.fastq.gz"
    tagged_fastq="$OUTPUT_DIR/${pre}.tagged.fastq.gz"
    fastq_screen --threads 16  --outdir $OUTPUT_DIR --force --tag --subset 100000 $pf 
    fastqSpliter.py --taggedFastq $tagged_fastq --prefix $pre --outDir $OUTPUT_DIR     # grep #FQST tag to the end and find the tag
done 




