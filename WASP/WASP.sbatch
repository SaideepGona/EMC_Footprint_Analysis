#!/bin/bash
###########################################
#SBATCH --job-name=WASP
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --output=logs/run_WASP_${PREFIX}.out
#SBATCH	--error=logs/run_WASP_${PREFIX}.err
###########################################

module load gcc/6.2
module load python/3.5.2
module load STAR/2.6.1a

# submit with qsub -v PREFIX=$NAME job.pbs (see submit_R_script_batch.sh from heard_mouse_data)
LOC=/project2/lbarreiro/users/katie/WASP

# STAR - initial mapping
STAR --runThreadN 5 --genomeDir /group/referenceFiles/Homo_sapiens/STARgenome/hg19_UCSC_50bp --readFilesIn ${LOC}/fastq/${PREFIX}.fastq --outFileNamePrefix ${LOC}/WASP_output/star_initial/$PREFIX --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat

# install required packages if not already and update path
# pip install --user pysam
# pip install --user tables
# export PATH=$HOME/.local:$PATH

# use find_intersecting_snps.py to identify reads that have mapping biases - takes a really long time 
python3 /apps/software/gcc-6.2.0/WASP/0.3.3/mapping/find_intersecting_snps.py --is_sorted --output_dir ${LOC}/WASP_output/intersecting_snps --snp_index ${LOC}/WASP_output/hdf5/snp_index.h5 --snp_tab ${LOC}/WASP_output/hdf5/snp_tab.h5 --haplotype ${LOC}/WASP_output/hdf5/haps.h5 --samples ${LOC}/WASP_inputs/sample_names.txt ${LOC}/WASP_output/star_initial/${PREFIX}Aligned.sortedByCoord.out.bam

# remapping the reads identified by find_intersecting_snps
STAR --runThreadN 5 --genomeDir /group/referenceFiles/Homo_sapiens/STARgenome/hg19_UCSC_50bp --readFilesIn ${LOC}/WASP_output/intersecting_snps/${PREFIX}Aligned.sortedByCoord.out.remap.fq.gz --outFileNamePrefix ${LOC}/WASP_output/star_remap/$PREFIX --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat

# filter out reads that don't remap to the same location as the original read
python3 /apps/software/gcc-6.2.0/WASP/0.3.3/mapping/filter_remapped_reads.py ${LOC}/WASP_output/intersecting_snps/${PREFIX}Aligned.sortedByCoord.out.to.remap.bam ${LOC}/WASP_output/star_remap/${PREFIX}Aligned.sortedByCoord.out.bam ${LOC}/WASP_output/filtered_remapped_reads/${PREFIX}.keep.bam

# merge the .keep.bam's and the .remap.keep.bam's to get a complete set of good reads
samtools merge ${LOC}/WASP_output/merged_bams/${PREFIX}.merge.bam ${LOC}/WASP_output/filtered_remapped_reads/${PREFIX}.keep.bam ${LOC}/WASP_output/intersecting_snps/${PREFIX}Aligned.sortedByCoord.out.keep.bam

samtools sort -o ${LOC}/WASP_output/merged_bams/${PREFIX}.merge.sort.bam ${LOC}/WASP_output/merged_bams/${PREFIX}.merge.bam

samtools index ${LOC}/WASP_output/merged_bams/${PREFIX}.merge.sort.bam

# filter duplicate reads
python3 /apps/software/gcc-6.2.0/WASP/0.3.3/mapping/rmdup.py ${LOC}/WASP_output/merged_bams/${PREFIX}.merge.sort.bam ${LOC}/WASP_output/duplicates_removed/${PREFIX}.nodup.bam	


