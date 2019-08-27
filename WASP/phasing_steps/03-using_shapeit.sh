#!/bin/bash
for chr in {1..22}; do

	IN_DIR=/project2/lbarreiro/users/katie/admixture/vcf_files/split_by_chr
	FILE=all_my_samples_merged_chr${chr}
	OUT_DIR=/project2/lbarreiro/users/katie/WASP

	/project2/lbarreiro/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit \
		--input-ped /project2/lbarreiro/users/katie/admixture/vcf_files/split_by_chr/all_my_samples_merged_chr${chr}.ped /project2/lbarreiro/users/katie/admixture/vcf_files/split_by_chr/all_my_samples_merged_chr${chr}.map \
		#--input-ped ${IN_DIR}/${FILE}.ped ${IN_DIR}/${FILE}.map \
		−−input−map genetic_map_GRCh37_chr${chr}.txt \
		--output-max /project2/lbarreiro/users/katie/WASP/phasing_outputs/all_my_samples_merged_chr${chr}.phased \
		--thread 8 \
		--output-log ${OUT_DIR}/phasing_outputs/${FILE}.phased
done
