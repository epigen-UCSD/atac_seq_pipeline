# activate conda environ
#source activate bds_atac_py3 

sample_prefix=$1
FASTQDIR="/projects/ps-epigen/seqdata/"
TMP_DIR="/oasis/tscc/scratch/$(whoami)/.tmp/"
mkdir -p $TMP_DIR

OUTPUT_DIR="/projects/ps-epigen/outputs/libQCs/${sample_prefix}/"
mkdir -p $OUTPUT_DIR


for i in `seq 1 2` 
do
    p="${sample_prefix}_R${i}.fastq"
    bz2file="${p}.bz2"
    echo "running decompress $bz2file ..."
    bzip2 -d -c $FASTQDIR$bz2file > $TMP_DIR$p

    echo "running fastqc $p ..."
    fastqc -t 16 -o $OUTPUT_DIR "${TMP_DIR}$p"

    fastq_screen --threads 16  --outdir $OUTPUT_DIR --force --tag --subset 100000  "${TMP_DIR}$p"
  
    tagged_fastq="$OUTPUT_DIR/${sample_prefix}_R${i}.tagged.fastq"
    [[ -f ${tagged_fastq}.gz ]] && rm ${tagged_fastq}.gz
    
    gzip -9 $tagged_fastq
    fastqSpliter.py --taggedFastq ${tagged_fastq}.gz --prefix ${sample_prefix}_R${i} --outDir $OUTPUT_DIR
    
done 

rm $TMP_DIR${sample_prefix}_*.fastq
