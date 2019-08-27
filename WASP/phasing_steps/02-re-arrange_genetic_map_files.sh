#!/bin/bash

for chr in {1..22}; do
    awk '{print $4, $2, $3}' genetic_map_GRCh37_chr$chr.txt > genetic_map_GRCh37_chr$chr.ordered.txt 
done