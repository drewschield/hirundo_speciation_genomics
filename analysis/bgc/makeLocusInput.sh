#!/bin/bash
pop=$1
while read region; do
	chrom=`echo "$region" | cut -f 1`;
	locus=`echo "$region" | cut -f 2`;
	echo formatting input data for $locus
	paste popmap.$pop.txt <(tail -n+2 ../hzar/input_background_hzar/data.$pop.$locus.txt | cut -f6) > ./input_locus/data.$pop.$locus.txt
done < background.$pop.snps.maf1.region.locus-table.txt

