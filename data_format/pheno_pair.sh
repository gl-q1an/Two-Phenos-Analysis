#!/bin/bash

# Generate Pheno pairs according to the pheno in filelist 
#   for subsequent batch analysis of genetic correlation.

set -e

filelist=filelist_template.txt
output_pairs=pheno_pairs.txt

ALL_PHENO=$(grep -v '^#' $filelist | cut -f 1 -d $'\t')
# MAIN_PHENO=$ALL_PHENO
MAIN_PHENO=PTSD


echo -e "# PHENO1\tPHENO2" > ${output_pairs}

for phe1 in $MAIN_PHENO;do
    if ! grep -q "^${phe1}[[:space:]]" "$filelist"; then
        echo "Warning: ${phe1} is not in the filelist!"
    fi
    for phe2 in $ALL_PHENO;do
        if [[ "${phe1}" != "${phe2}" ]]; then
            if ! grep -q -P "^$phe2\t$phe1$" "${output_pairs}"; then
                echo -e "${phe1}\t${phe2}" >> ${output_pairs}
            fi
        fi
    done
done