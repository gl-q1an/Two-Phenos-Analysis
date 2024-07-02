#!/bin/bash

set -e

export PLEIOFDR_DIR=/home/user/genotools/py_cov
export PLEIOFDR_REF=/home/user/reference/pleiofdr
export PYCOV_REF=/home/user/reference/py_cov
export filelist=../data_format/filelist_template.txt
export output_config_dir=./pleioFDR/config
export output_dir=./pleioFDR/result
export activate=/home/user/miniconda3/bin/activate
export env=py2

# Config File Variable
export template_cfg=${PLEIOFDR_DIR}/config_template.txt

export ref_file=${PLEIOFDR_REF}/ref9545380_1kgPhase3eur_LDr2p1.mat
export trait_folder=./pleioFDR/format
export format_mat="@.mat"
export refinfo=${PLEIOFDR_REF}/9545380.ref
export randprune_n=200
export exclude_chr_pos="6 25119106 33854733" 

cFDR_run(){
    local pheno1=$1
    local pheno2=$3

    local out_config_dir_cfdr=${output_config_dir}/${pheno1}_${pheno2}
    local out_result_dir_cfdr=${output_dir}/${pheno1}_${pheno2}

    mkdir -p ${out_config_dir_cfdr}
    mkdir -p ${out_result_dir_cfdr}

    ##### The First Part : Write the config #####
    # confFDR and cond1FDR
    local trait_name1=${pheno1}
    local trait_file1=$(echo ${format_mat} | sed "s/@/${pheno1}/")
    local trait_name2=${pheno2}
    local trait_file2=$(echo ${format_mat} | sed "s/@/${pheno2}/")

    stat_type=conjfdr
    fdr_thresh=0.05
    out_dir=${out_result_dir_cfdr}/conjFDR
    >${out_config_dir_cfdr}/conj_config.txt
    while read line;do
        eval "echo \"$line\"" >> ${out_config_dir_cfdr}/conj_config.txt
    done < ${template_cfg}

    stat_type=condfdr
    fdr_thresh=0.01
    out_dir=${out_result_dir_cfdr}/cond1FDR
    >${out_config_dir_cfdr}/cond1_config.txt
    while read line;do
        eval "echo \"$line\"" >> ${out_config_dir_cfdr}/cond1_config.txt
    done < ${template_cfg}

    # cond2FDR
    local trait_name1=${pheno2}
    local trait_file1=$(echo ${format_mat} | sed "s/@/${pheno2}/")   
    local trait_name2=${pheno1}
    local trait_file2=$(echo ${format_mat} | sed "s/@/${pheno1}/")

    stat_type=condfdr
    fdr_thresh=0.01
    out_dir=${out_result_dir_cfdr}/cond2FDR
    >${out_config_dir_cfdr}/cond2_config.txt
    while read line;do
        eval "echo \"$line\"" >> ${out_config_dir_cfdr}/cond2_config.txt
    done < ${template_cfg}
    ##### The First Part End!#####

    ##### RUN #####
    conjconfig=${out_config_dir_cfdr}/conj_config.txt
    cond1config=${out_config_dir_cfdr}/cond1_config.txt
    cond2config=${out_config_dir_cfdr}/cond2_config.txt
    matlab -nodisplay -nosplash -r "config='${conjconfig}'; run('${pleiofdr_dir}/runme'); exit"
    matlab -nodisplay -nosplash -r "config='${cond1config}'; run('${pleiofdr_dir}/runme'); exit"
    matlab -nodisplay -nosplash -r "config='${cond2config}'; run('${pleiofdr_dir}/runme'); exit"
}

export -f cFDR_run

grep -v '^#' $pairlist |\
    xargs -n 2 bash -c 'cFDR_run "$0" "$1"'

# if you want run it parallelly
if false; then
    export jobs=10
    grep -v '^#' $pairlist |\
        parallel --jobs ${jobs} --colsep "\t" \
        'cFDR_run {1} {2}'
fi