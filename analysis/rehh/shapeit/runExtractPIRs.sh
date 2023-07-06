chromlist=$1
for chrom in `cat chrom.list`; do
	./extractPIRs/extractPIRs --bam bamlist.$chrom --vcf ../input/hirundo_rustica.parental.$chrom.snps.miss05.maf01.vcf --out ./pirs/pirs.$chrom --base-quality 20 --read-quality 20
done

