#!/bin/bash

set -e

export LDSC_DIR=/home/user/genotools/ldsc
export LDSC_REF_DIR=/home/user/reference/ldsc
export signed="Z,0" # If z does not exist in format, you can use BETA,0 or OR,1
export ldsc_sum_dir="/home/user/data/ldsc"
export ldsc_sum="${ldsc_sum_dir}/@.sumstats.gz"
export pairlist=../data_format/pheno_pairs.txt
export output_dir=./ldsc
export activate=/home/user/miniconda3/bin/activate

mkdir -p ${output_dir}

ldsc_corr(){
    local input_phe1=$1
    local input_phe2=$2
    local ldsc_ref=$3
    local outdir=$4
    mkdir -p ${output_dir}
    source ${activate} ldsc
    ${LDSC_DIR}/ldsc.py \
        --rg ${input_phe1},${input_phe2} \
        --ref-ld-chr ${ldsc_ref}/ \
        --w-ld-chr ${ldsc_ref}/ \
        --out ${outdir}/${input_phe1}_${input_phe2}_ldsc
}

export -f ldsc_corr

grep -v '^#' $pairlist |\
    xargs -n 2 bash -c 'ldsc_corr "${ldsc_sum//@/$0}" "${ldsc_sum//@/$1}" "${LDSC_REF_DIR}/eur_w_ld_chr" "${output_dir}"' 

# if you want run it parallelly
if false; then
    export jobs=10
    grep -v '^#' $pairlist |\
        parallel --jobs ${jobs} --colsep "\t" \
        ldsc_format "$(echo ${format_sum} | sed "s/@/{1}/")" "$(echo ${format_sum} | sed "s/@/{2}/")" "${LDSC_REF_DIR}/eur_w_ld_chr" "${output_dir}"
fi