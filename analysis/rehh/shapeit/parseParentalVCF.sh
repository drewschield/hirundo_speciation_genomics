table=$1
while read i; do 
	scaff=`echo "$i" | cut -f 1`
	chrom=`echo "$i" | cut -f 2`
	echo "parsing parental genotypes for $chrom"
	bcftools view --threads 16 -S popmap.all -t $scaff -q 0.01:minor -m2 -M2 -U -v snps -i 'F_MISSING<0.5' -O z -o ./input/hirundo_rustica.parental.$chrom.snps.miss05.maf01.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.ingroup.vcf.gz
done < $table

