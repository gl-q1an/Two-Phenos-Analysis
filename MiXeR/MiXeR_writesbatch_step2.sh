#!/bin/bash

set -e

export MIXER_DIR=/home/user/genotools/mixer
export MIXER_REF_DIR=/home/user/reference/mixer
export MIXER_FORMAT_DIR=./MiXeR/format
export fit1_dir=./MiXeR/fit1
export pairlist=../data_format/pheno_pairs.txt
export activate=/home/user/miniconda3/bin/activate
export env=MiXeR

# WRITE THE FIT1 FILE
export output_dir=./MiXeR/sbatchfile/ftc2
export ftc2_template=sbatch_template_step2.txt

export ftc2_workdir=./MiXeR/ftc2

mkdir -p ${output_dir}

write_sbatch_ftc2(){
    local pheno1=$1
    local pheno2=$2
    local sbatchfile=${output_dir}/ftc2_${pheno1}_${pheno2}.sbatch
    cp ${ftc2_template} ${sbatchfile}
    sed -i "s|<PHENO1_TO_SWAP>|${pheno1}|g" ${sbatchfile}
    sed -i "s|<PHENO2_TO_SWAP>|${pheno2}|g" ${sbatchfile}
    sed -i "s|<MIXERDIR_TO_SWAP>|${MIXER_DIR}|g" ${sbatchfile}
    sed -i "s|<WORKDIR_TO_SWAP>|${ftc2_workdir}|g" ${sbatchfile}
    sed -i "s|<FIT1DIR_TO_SWAP>|${fit1_dir}|g" ${sbatchfile}
    sed -i "s|<MIXER_REF_DIR_TO_SWAP>|${MIXER_REF_DIR}|g" ${sbatchfile}
    sed -i "s|<FORMAT_DIR_TO_SWAP>|${MIXER_FORMAT_DIR}|g" ${sbatchfile}
    sed -i "s|<ACTIVATE_DIR_TO_SWAP>|${activate}|g" ${sbatchfile}
    sed -i "s|<ENV_DIR_TO_SWAP>|${env}|g" ${sbatchfile}
}

grep -v '^#' $pairlist |\
    xargs -n 2 bash -c 'write_sbatch_ftc2 "$0" "$1"'