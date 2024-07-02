#!/bin/bash

set -e

export MIXER_DIR=/home/user/genotools/mixer
export MIXER_REF_DIR=/home/user/reference/mixer
export MIXER_FORMAT_DIR=./MiXeR/format
export filelist=../data_format/filelist_template.txt
export activate=/home/user/miniconda3/bin/activate
export env=MiXeR

# WRITE THE FIT1 FILE
export output_dir=./MiXeR/sbatchfile/fit1
export fit1_template=sbatch_template_step1.txt

export fit1_workdir=./MiXeR/fit1

mkdir -p ${output_dir}

write_sbatch_fit1(){
    local pheno=$1
    local sbatchfile=${output_dir}/fit1_${pheno}.sbatch
    cp ${fit1_template} ${sbatchfile}
    sed -i "s|<PHENO_TO_SWAP>|${pheno}|g" ${sbatchfile}
    sed -i "s|<MIXERDIR_TO_SWAP>|${MIXER_DIR}|g" ${sbatchfile}
    sed -i "s|<WORKDIR_TO_SWAP>|${fit1_workdir}|g" ${sbatchfile}
    sed -i "s|<MIXER_REF_DIR_TO_SWAP>|${MIXER_REF_DIR}|g" ${sbatchfile}
    sed -i "s|<FORMAT_DIR_TO_SWAP>|${MIXER_FORMAT_DIR}|g" ${sbatchfile}
    sed -i "s|<ACTIVATE_DIR_TO_SWAP>|${activate}|g" ${sbatchfile}
    sed -i "s|<ENV_DIR_TO_SWAP>|${env}|g" ${sbatchfile}
}

grep -v '^#' $filelist |\
    xargs -n 2 bash -c 'write_sbatch_fit1 "$0"'