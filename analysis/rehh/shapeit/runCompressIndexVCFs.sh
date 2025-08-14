chromlist=$1
for chrom in `cat chrom.list`; do
	bgzip ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased.vcf
	tabix -p vcf ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased.vcf.gz
done
