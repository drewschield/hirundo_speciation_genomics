chromlist=$1
for chrom in `cat chrom.list`; do
	shapeit2/bin/shapeit -convert --input-haps ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased --output-vcf ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased.vcf
done

