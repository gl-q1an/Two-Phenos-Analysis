#!/bin/bash

set -e

export LAVA_RSCRIPT="$(dirname "$0")"/func_lava.R
export LAVA_CONFIG=lava_config_template.txt
export LAVA_REF_DIR=/home/user/reference/lava
export signed="Z,0" # If z does not exist in format, you can use BETA,0 or OR,1
export format_sum_dir="/home/user/data/format"
export format_sum="${format_sum_dir}/@_formated.txt.gz"
export pairlist=../data_format/pheno_pairs.txt
export output_dir=./lava
export activate=/home/user/miniconda3/bin/activate
export env=R4.2
export jobs=10

mkdir -p ${output_dir}

lava_para(){
    local config=$1
    local phe1=$2
    local phe2=$3
    local outfile=$4
    local refdir=$5
    local max=$6
    source ${act} ${env}
    echo "${phe1} and ${phe2} begin at $(date)"
    mkdir -p $(dirname ${outfile})
    Rscript ${script_dir}/func_lava.R \
    --lavaconfig "${config}" \
    --phe1 ${phe1} \
    --phe2 ${phe2} \
    --outfile ${outfile} \
    --refdir ${refdir} \
    --max ${max}
    echo "${phe1} and ${phe2} end at $(date)"
}

grep -v '^#' $pairlist |\
    xargs -n 2 bash -c 'lava_para "${LAVA_CONFIG}" "$0" "$1" "${output_dir}/lava_${0}_${1}/lava_${0}_${1}" "${LAVA_REF_DIR}" "${jobs}"' 

# Because Lava comes with parallel processing, so we don't provide parallel mode here.