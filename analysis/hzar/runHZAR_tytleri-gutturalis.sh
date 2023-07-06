while read region; do
	chrom=`echo "$region" | cut -f 1`;
	locus=`echo "$region" | cut -f 2`;
	echo performing geographic cline analysis on $locus
	Rscript hzar_background_tytleri-gutturalis.R $locus
done < background.tytleri-gutturalis.snps.maf1.region.locus-table.txt
