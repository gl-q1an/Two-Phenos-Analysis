#!/bin/bash

set -e

export PYCOV_DIR=/home/user/genotools/py_cov
export format_sum_dir="/home/user/data/format"
export format_sum="${format_sum_dir}/@_formated.txt.gz"
export filelist=../data_format/filelist_template.txt
export output_dir=./MiXeR/format
export activate=/home/user/miniconda3/bin/activate
export env=MiXeR

mixer_format(){
    local inputsum=$1
    local outfile_pre=$2

    source ${activate} ${env}
    mkdir -p $(dirname ${outfile_pre})

    python ${PYCOV_DIR}/sumstats.py csv \
        --sumstats ${inputsum} \
        --force --auto --head 5 \
        --ignore BETA \
        --out ${outfile_pre}.csv

    python ${PYCOV_DIR}/sumstats.py zscore \
        --sumstats ${outfile_pre}.csv |\
    python ${PYCOV_DIR}/sumstats.py qc \
        --exclude-ranges 6:26000000-34000000 \
        --out ${outfile_pre}_noMHC.csv --force

    rm -f ${outfile_pre}.csv

    gzip ${outfile_pre}_noMHC.csv
    echo "$(basename ${outfile_pre}) is formated!!"
}

export -f mixer_format

grep -v '^#' $filelist |\
    xargs -n 2 bash -c 'mixer_format "${format_sum//@/$0}" "${output_dir}/${0}/${0}"' 

# if you want run it parallelly
if false; then
    export jobs=10
    grep -v '^#' $filelist | cut -f 1 -d $'\t' |\
        parallel --jobs ${jobs} mixer_format "$(echo ${format_sum} | sed "s/@/{}/")" "${output_dir}/{}/{}"
fi