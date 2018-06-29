#!/bin/bash
#module load fastqc

sample_prefix=$1
FASTQDIR="/projects/ps-epigen/seqdata/"
OUTPUT_DIR="/projects/ps-epigen/outputs/libQCs/${sample_prefix}"

mkdir -p $OUTPUT_DIR

for i in `seq 1 2`
do
    pre="${sample_prefix}_R${i}.trim"
    
    [[ ! -f "${FASTQDIR}${pre}.fastq.gz" ]] &&  pre=${pre/.trim/} 

    tagged_fastq="$OUTPUT_DIR/${pre}.tagged.fastq.gz"
    pf="${FASTQDIR}${pre}.fastq.gz"
    
    echo "running fastqc $pf ..."
    fastqc  -t 8 -o $OUTPUT_DIR $pf

    echo "running fastscreen $p ..."    
    fastq_screen --threads 8  --outdir $OUTPUT_DIR --force --tag --subset 100000 $pf 


    # grep #FQST tag to the end and find the tag
    fastqSpliter.py --taggedFastq $tagged_fastq --prefix $pre --outDir $OUTPUT_DIR
done 




