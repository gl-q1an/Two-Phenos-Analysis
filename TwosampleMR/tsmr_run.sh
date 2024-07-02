#!/bin/bash

set -e

export TSMR_RSCRIPT="$(dirname "$0")"/func_tsmr.R
export format_sum_dir="/home/user/data/format"
export format_sum="${format_sum_dir}/@_formated.txt.gz"
export pairlist=../data_format/pheno_pairs.txt
export PLINK_PATH=/home/user/genotools/plink/plink
export REF_1KG=/home/user/reference/magma/g1000_eur/g1000_eur
export output_dir=./tsmr/result
export p_thresh=5e-8
export activate=/home/user/miniconda3/bin/activate
export env=R4.2

mkdir -p ${output_dir}

tsmr_run(){
    local pheno1=$1
    local pheno2=$2
    local pheno1sum=$(echo ${format_sum} | sed "s/@/${1}/")
    local pheno2sum=$(echo ${format_sum} | sed "s/@/${2}/")
    source ${activate} ${env}
    Rscript ${TSMR_RSCRIPT} \
        -P ${PLINK_PATH} \
        -R ${REF_1KG} \
        -o ${output_dir} \
        -x ${pheno1} \
        -f ${pheno1sum} \
        -y ${pheno2} \
        -F ${pheno2sum} \
        -p 5e-8 \
        -t YourLinkRToken
}

grep -v '^#' $pairlist |\
    xargs -n 2 bash -c 'tsmr_run "$0" "$1"'

# LDlinkR can only use one API to search at a time, 
#    so this step cannot be performed in parallel.