#!/bin/bash
# setQC1: get the multiQC reports for the input libs

usage() { echo "Usage: $0 [-l <ints>] [-s <int>]" 1>&2; exit 1; }

while getopts "l:s:" o; do
    case $o in
        l) l+=("${OPTARG}");;
        s) s="${OPTARG}";;
        *) usage;;
    esac
done

shift $((OPTIND-1))

#if [ -z "${l}" ] || [ -z "${s}" ]; then
#    usage
#fi

echo "input lib number: ${l[@]}"
echo "output set number = ${s}"

# load the multiQC environment
source activate bds_atac_py3


# take the input as array
#libs_array=( "$(@)" )
libs_array=(${l})
set_no=$s
libQC_dir="/projects/ps-epigen/outputs/"
out_dir="/projects/ps-epigen/outputs/setQCs/Set_${set_no}/"
mkdir -p $out_dir

echo "${libs_array[@]/#/JYH_}"> ${out_dir}including_libs.txt

# paste the command 
cmd="multiqc -k tsv -f -p ${libs_array[@]/#/${libQC_dir}libQCs/JYH_}  ${libs_array[@]/#/${libQC_dir}peaks/JYH_} -o $out_dir  "

echo "running command:"
echo $cmd 
eval $cmd 










