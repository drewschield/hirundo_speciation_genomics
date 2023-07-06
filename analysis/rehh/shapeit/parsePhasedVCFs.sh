chromlist=$1
for chrom in `cat chrom.list`; do
	echo parsing $chrom
	bcftools view --threads 16 -S popmap.rustica -O z -o ./results/hirundo_rustica.parental.rustica.$chrom.snps.miss05.maf01.phased.vcf.gz ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased.vcf.gz
	bcftools view --threads 16 -S popmap.tytleri -O z -o ./results/hirundo_rustica.parental.tytleri.$chrom.snps.miss05.maf01.phased.vcf.gz ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased.vcf.gz
	bcftools view --threads 16 -S popmap.gutturalis -O z -o ./results/hirundo_rustica.parental.gutturalis.$chrom.snps.miss05.maf01.phased.vcf.gz ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased.vcf.gz
done
