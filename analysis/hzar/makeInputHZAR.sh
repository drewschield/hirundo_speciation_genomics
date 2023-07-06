#!/bin/bash
pop=$1
counter=1
while read region; do
	chrom=`echo "$region" | cut -f 1`;
	start=`echo "$region" | cut -f 2`;
	end=`echo "$region" | cut -f 3`;
	bed=`echo "$chrom:$start-$end"`
	echo -e "$bed\tlocus$counter" >> background.$pop.snps.maf1.region.locus-table.txt
	paste data.$pop.txt <(cut -f2 ./results_introgress/$pop.locus.$bed.txt) > ./input_background_hzar/data.$pop.locus$counter.txt
	counter=$((counter+1))
done < background.$pop.snps.maf1.region.txt
