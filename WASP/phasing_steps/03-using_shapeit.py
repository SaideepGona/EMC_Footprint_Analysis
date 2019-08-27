import os

command = [
    "/project2/lbarreiro/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit",
	"--input-ped", "place2", "place3",
	"--input-map", "place5",
	"--output-max", "place7",
	"--thread", "1",
	"--output-log", "place11",
]

for chr in range(1,23):

	IN_DIR="/project2/lbarreiro/users/katie/admixture/vcf_files/split_by_chr_with_mind_flag/"
	FILE="all_my_samples_merged_chr" + str(chr)
	OUT_DIR="/project2/lbarreiro/users/katie/WASP"

	command_cur = command[:]

	command_cur[2] = IN_DIR + FILE + ".ped"
	command_cur[3] = IN_DIR + FILE + ".map"
	command_cur[5] = "genetic_map_GRCh37_chr" + str(chr) + ".ordered.txt"
	command_cur[7] = OUT_DIR + "/phasing_outputs/" + FILE + ".phased"
	command_cur[11] = OUT_DIR + "/phasing_outputs/" + FILE + ".phased.log"

	command_cur = " ".join(command_cur)

	print(command_cur)

	os.system(command_cur)
	


# 	/project2/lbarreiro/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit \
# 		--input-ped /project2/lbarreiro/users/katie/admixture/vcf_files/split_by_chr/all_my_samples_merged_chr${chr}.ped /project2/lbarreiro/users/katie/admixture/vcf_files/split_by_chr/all_my_samples_merged_chr${chr}.map \
# 		#--input-ped ${IN_DIR}/${FILE}.ped ${IN_DIR}/${FILE}.map \
# 		−−input−map genetic_map_GRCh37_chr${chr}.txt \
# 		--output-max /project2/lbarreiro/users/katie/WASP/phasing_outputs/all_my_samples_merged_chr${chr}.phased \
# 		--thread 8 \
# 		--output-log ${OUT_DIR}/phasing_outputs/${FILE}.phased
# done
