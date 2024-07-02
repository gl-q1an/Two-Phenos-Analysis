#!/bin/bash

set -e

export LDSC_DIR=/home/user/genotools/ldsc
export LDSC_REF_DIR=/home/user/reference/ldsc
export signed="Z,0" # If z does not exist in format, you can use BETA,0 or OR,1
export format_sum_dir="/home/user/data/format"
export format_sum="${format_sum_dir}/@_formated.txt.gz"
export filelist=../data_format/filelist_template.txt
export output_dir=./ldsc
export activate=/home/user/miniconda3/bin/activate

mkdir -p ${output_dir}

ldsc_format(){
    local input=$1
    local output=$2
    local ref=$3
    source ${activate} ldsc
    if [ -n "$ref" ]; then
        ${LDSC_DIR}/munge_sumstats.py \
             --sumstats "${input}" \
             --merge-alleles "${ref}" \
             --signed-sumstats ${signed} \
             --out "${output}"
    else
        ${LDSC_DIR}/munge_sumstats.py \
             --sumstats "${input}" \
             --merge-alleles "${ref}" \
             --signed-sumstats ${signed} \
             --out "${output}"
    fi
}

export -f ldsc_format

grep -v '^#' $filelist |\
    xargs -n 2 bash -c 'ldsc_format "${format_sum//@/$0}" "${output_dir}/${0}"' 

# if you want run it parallelly
if false; then
    export jobs=10
    grep -v '^#' $filelist | cut -f 1 -d $'\t' |\
        parallel --jobs ${jobs} ldsc_format "$(echo ${format_sum} | sed "s/@/{}/")" "${output_dir}/{}"
fi