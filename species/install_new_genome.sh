#!/bin/bash

GENOME=$1 #rn6
target_dir=$2 # ~/data/GENOME
target_dir=${target_dir}/${GENOME} # ~/data/GENOME
SPECIES_FILE="/home/zhc268/data/software/atac_dnase_pipelines/species/epigen.conf"

ucsc_base_url="rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/"
genome_url=$ucsc_base_url${GENOME}.fa.gz


mkdir -p $target_dir 

# 1. download
echo "downloading seqences from ucsc" 
rsync -avzP $genome_url $target_dir



# 2. seq folder/ extract fasta per chromosome
# source activate bds_atac
cd $target_dir
echo "Extracting/processing data files..."

REF_FA_PREFIX=$GENOME".fa"
gzip -d -f -c ${REF_FA_PREFIX}.gz > ${REF_FA_PREFIX}

mkdir -p seq
cd seq
#rm -f ${REF_FA_PREFIX}
ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
faidx -x ${REF_FA_PREFIX}
cp --remove-destination *.fai ../

## determine gensz
CHRSZ=$GENOME.chrom.sizes
cut -f1,2 ${REF_FA_PREFIX}.fai > ../$CHRSZ

cd $target_dir
GENSZ=$(cat $CHRSZ | awk '{sum+=$2} END{print sum}')

## bowtie2_index
mkdir bowtie2_index
if [ ! -f ./bowtie2_index/ ${REF_FA_PREFIX}.rev.1.bt2 ]; then
    bowtie2-build ${REF_FA_PREFIX} ${REF_FA_PREFIX}
    mv *.bt2 ./bowtie2_index/
    ln -s  ${REF_FA_PREFIX}.fa ./bowtie2_index/
fi



# 3. add to the species file
if [[ $UMAP != "" ]]; then UMAP_PATH="${DATA_DIR}/$GENOME/$(basename $UMAP .tgz)"; fi
if [[ $BLACKLIST != "" ]]; then BLACKLIST_PATH="${DATA_DIR}/$GENOME/$(basename $BLACKLIST)"; fi
if [[ $TSS_ENRICH != "" ]]; then TSS_ENRICH_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $TSS_ENRICH)"; fi
if [[ $DNASE != "" ]]; then DNASE_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $DNASE)"; fi
if [[ $PROM != "" ]]; then PROM_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $PROM)"; fi
if [[ $ENH != "" ]]; then ENH_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $ENH)"; fi
if [[ $REG2MAP != "" ]]; then REG2MAP_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $REG2MAP)"; fi
if [[ $REG2MAP_BED != "" ]]; then REG2MAP_BED_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $REG2MAP_BED)"; fi
if [[ $ROADMAP_META != "" ]]; then ROADMAP_META_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $ROADMAP_META)"; fi

echo -e "[$GENOME] # installed by install_genome_data.sh" >> ${SPECIES_FILE}
echo -e "chrsz\t= ${DATA_DIR}/$GENOME/$(basename $CHRSZ)" >> ${SPECIES_FILE}
echo -e "seq\t= ${DATA_DIR}/$GENOME/seq" >> ${SPECIES_FILE}
echo -e "gensz\t= $GENSZ" >> ${SPECIES_FILE}
if [[ $UMAP != "" ]]; then echo -e "umap\t= ${UMAP_PATH}" >> ${SPECIES_FILE}; fi
if [ ${BUILD_BWT2_IDX} == 1 ]; then
    echo -e "bwt2_idx\t= ${DATA_DIR}/$GENOME/bowtie2_index/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
fi

if [ ${BUILD_BWA_IDX} == 1 ]; then
    echo -e "bwa_idx\t= ${DATA_DIR}/$GENOME/bwa_index/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
fi

echo -e "ref_fa\t= ${DATA_DIR}/$GENOME/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
if [[ ${BLACKLIST_PATH} != "" ]]; then echo -e "blacklist\t= ${BLACKLIST_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${SPECIES_BROWSER} != "" ]]; then echo -e "species_browser\t= ${SPECIES_BROWSER}" >> ${SPECIES_FILE}; fi
if [[ ${TSS_ENRICH_PATH} != "" ]]; then echo -e "# data for ATAQC" >> ${SPECIES_FILE}; fi
if [[ ${TSS_ENRICH_PATH} != "" ]]; then echo -e "tss_enrich\t= ${TSS_ENRICH_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${DNASE_PATH} != "" ]]; then echo -e "dnase\t= ${DNASE_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${PROM_PATH} != "" ]]; then echo -e "prom\t= ${PROM_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${ENH_PATH} != "" ]]; then echo -e "enh\t= ${ENH_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${REG2MAP_PATH} != "" ]]; then echo -e "reg2map\t= ${REG2MAP_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${REG2MAP_BED_PATH} != "" ]]; then echo -e "reg2map_bed\t= ${REG2MAP_BED_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${ROADMAP_META_PATH} != "" ]]; then echo -e "roadmap_meta\t= ${ROADMAP_META_PATH}" >> ${SPECIES_FILE}; fi
if [[ ${EXTRA_LINE} != "" ]]; then echo -e "${EXTRA_LINE}" >> ${SPECIES_FILE}; fi
echo "" >> ${SPECIES_FILE}


grep -Ev "(random)|(chrUn)" rn6_tss.bed.txt | awk '{print $1}' | sort | uniq
grep -Ev "(random)|(chrUn)" rn6_tss.bed.txt | wc -l 
                                                     




                                                               
