#!/bin/bash

set -e

export PYCOV_DIR=/home/user/genotools/py_cov
export PYCOV_REF=/home/user/reference/py_cov
export format_sum_dir="/home/user/data/format"
export format_sum="${format_sum_dir}/@_formated.txt.gz"
export filelist=../data_format/filelist_template.txt
export output_dir=./pleioFDR/format
export activate=/home/user/miniconda3/bin/activate
export env=py2

cFDR_format(){
    local inputsum=$1
    local outfile=$2
    local ref=$3

    mkdir -p $(dirname "${outfile}")
    source ${activate} ${env}
    python ${PYCOV_DIR}/sumstats.py csv \
    --auto --force \
    --sumstats ${inputsum} \
    --ignore BETA \
    --out ${outfile}.csv

    python ${PYCOV_DIR}/sumstats.py zscore \
    --sumstats ${outfile}.csv \
    --out ${outfile}.zscore \
    --force

    python ${PYCOV_DIR}/sumstats.py mat \
    --sumstats ${outfile}.zscore \
    --ref ${ref} \
    --out ${outfile}.mat \
    --force
}

export -f cFDR_format

grep -v '^#' $filelist |\
    xargs -n 2 bash -c 'cFDR_format "${format_sum//@/$0}" "${output_dir}/${0}/${0}" "${PYCOV_REF}/9545380.ref"'

# if you want run it parallelly
if false; then
    export jobs=10
    grep -v '^#' $filelist | cut -f 1 -d $'\t' |\
        parallel --jobs ${jobs} cFDR_format "$(echo ${format_sum} | sed "s/@/{}/")" "${output_dir}/{}/{}" "${PYCOV_REF}/9545380.ref"
fi