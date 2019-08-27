#!/bin/bash

LOC=/project2/lbarreiro/users/katie/WASP_inputs/prefix_list.txt

while read LINE; do
	sbatch -v PREFIX=$LINE WASP.sbatch
done < $LOC

