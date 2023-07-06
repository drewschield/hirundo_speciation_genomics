chromlist=$1
for chrom in `cat chrom.list`; do
	shapeit2/bin/shapeit -assemble --thread 8 --input-vcf ../input/hirundo_rustica.parental.$chrom.snps.miss05.maf01.vcf --input-pir ./pirs/pirs.$chrom --states 1000 --burn 200 --prune 210 --main 2000 --force -O ./results/hirundo_rustica.parental.$chrom.snps.miss05.maf01.phased
done

