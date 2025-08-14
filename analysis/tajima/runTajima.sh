for pop in rustica tytleri gutturalis; do
	echo estimating Tajima D in $pop
	vk tajima 100000 100000 ./vcf/hirundo_rustica.parental.$pop.snps.vcf.gz > ./results/tajima.$pop.100kb.txt
	vk tajima 100000 100000 ./vcf/hirundo_rustica.parental.$pop.snps.chrZ.vcf.gz | tail -n+2 >> ./results/tajima.$pop.100kb.txt
	vk tajima 10000 10000 ./vcf/hirundo_rustica.parental.$pop.snps.vcf.gz > ./results/tajima.$pop.10kb.txt
	vk tajima 10000 10000 ./vcf/hirundo_rustica.parental.$pop.snps.chrZ.vcf.gz | tail -n+2 >> ./results/tajima.$pop.10kb.txt
done
