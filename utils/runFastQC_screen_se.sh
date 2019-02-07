# activate conda environ
#source activate bds_atac_py3 

sample_prefix=$1
FASTQDIR="/projects/ps-epigen/seqdata/"
OUTPUT_DIR="/projects/ps-epigen/outputs/libQCs/${sample_prefix}"
WORKDIR="/oasis/tscc/scratch/$(whoami)/outputs/"

mkdir -p $OUTPUT_DIR

# do fastqc on raw fastq files
pre="${sample_prefix}"
[[ ! -f "${FASTQDIR}${pre}.fastq.gz" ]] &&  pre=${pre}.trim
pf="${FASTQDIR}${pre}.fastq.gz"
echo "running fastqc $p ..."
fastqc -t 16 -o $OUTPUT_DIR $pf

# do fastq_screen on trimmed fastq files (perfer)
pre="${sample_prefix}.trim"
pf=$(find ${WORKDIR}$sample_prefix -name "${pre}.fastq.gz")
[[ -z $pf ]] &&  pre=${pre/.trim/} && pf="${FASTQDIR}${pre}.fastq.gz"
tagged_fastq="$OUTPUT_DIR/${pre}.tagged.fastq.gz"
echo "running fastscreen $pf ..."
fastq_screen --threads 16  --outdir $OUTPUT_DIR --force --tag --subset 100000 $pf
fastqSpliter.py --taggedFastq $tagged_fastq --prefix $pre --outDir $OUTPUT_DIR     # grep #FQST tag to the end and find the tag

    

