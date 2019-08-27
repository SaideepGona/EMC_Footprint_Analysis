#!/bin/bash

module load plink
FILE=/project2/lbarreiro/users/katie/admixture/vcf_files/all_my_samples_merged

for chr in {1..22}; do
     plink --bfile ${FILE} \
	   --chr $chr \
	   --allow-extra-chr \
	   --mind 0.1 \
	   --geno 0.01 \
	   --recode \
	   --out ${FILE}_chr$chr ;
done     
