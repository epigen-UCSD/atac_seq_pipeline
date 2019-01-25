# activate conda environ
#source activate bds_atac_py3 

sample_prefix=$1
FASTQDIR="/projects/ps-epigen/seqdata/"
OUTPUT_DIR="/projects/ps-epigen/outputs/libQCs/${sample_prefix}"
WORKDIR="/oasis/tscc/scratch/$(whoami)/outputs/"

mkdir -p $OUTPUT_DIR

#find $OUTPUT_DIR -name "*screen*" -delete

p="${sample_prefix}.trim.fastq.gz"
[[ ! -f "${FASTQDIR}${p}" ]] &&  p=${p/.trim/}

echo "running fastqc $p ..."
fastqc -t 16 -o $OUTPUT_DIR "${FASTQDIR}$p"
fastq_screen --threads 16  --outdir $OUTPUT_DIR --force --tag --subset 100000 "${FASTQDIR}$p"
  
tagged_fastq="$OUTPUT_DIR/${sample_prefix}.trim.tagged.fastq"
[[ -f ${tagged_fastq}.gz ]] && rm ${tagged_fastq}.gz
gzip -9 $tagged_fastq
fastqSpliter.py --taggedFastq ${tagged_fastq}.gz --prefix ${sample_prefix} --outDir $OUTPUT_DIR
    

