pop=$1
while read region; do
	chrom=`echo "$region" | cut -f 1`;
	snp=`echo "$region" | cut -f 2`;
	tmp=`grep $chrom all.$pop.snps.maf1.txt | grep -w -B49 -A50 $snp`
	start=`echo $tmp | cut -d' ' -f2`
	end=`echo $tmp | rev | cut -d' ' -f 1 | rev`
	echo "$chrom\t$start\t$end" >> background.$pop.snps.maf1.region.txt
done < background.$pop.snps.maf1.txt

