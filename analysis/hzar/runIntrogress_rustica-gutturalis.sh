while read region; do
	chrom=`echo "$region" | cut -f 1`;
	start=`echo "$region" | cut -f 2`;
	end=`echo "$region" | cut -f 3`;
	bed=`echo "$chrom:$start-$end"`
	echo "estimating hybrid indices for ${bed}"
	Rscript introgress_rustica-gutturalis.R $bed 
done < background.rustica-gutturalis.snps.maf1.region.txt
