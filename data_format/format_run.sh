#!/bin/bash

set -e

export form=SNP,CHR,BP,A1,A2,Frq,P,BETA,SE,N
export filelist=filelist_template.txt
export output_dir=

mkdir -p ${output_dir}

format_summary(){
    local summary=$1
    local output=$2 
    python "$(dirname "$0")"/format_summary.py \
      --summary "${summary}" \
      --output "${output}" \
      --outform "${form}" 
}

export -f format_summary

grep -v '^#' $filelist |\
    xargs -n 2 bash -c 'funmat_summary "$1" "${output_dir}/${0}_formated.txt.gz"' 

# if you want run it parallelly
if false; then
    export jobs=10

    grep -v '^#' $filelist |\
        parallel --jobs ${jobs} --colsep "\t" format_summary {2} "${output_dir}/{1}_formated.txt.gz"
fi