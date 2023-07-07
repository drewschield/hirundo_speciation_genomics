for pop in rustica tytleri gutturalis rustica-tytleri rustica-gutturalis tytleri-gutturalis; do
	vcftools --gzvcf /data3/hirundo/analysis/ld/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.phased.vcf.gz --keep popmap.rand.$pop --bed Hirundo_rustica_bHirRus1.final.window_100kb.auto.rand.bed --maf 0.05 --hap-r2 --ld-window-bp 25000 --out ./r2/r2.$pop.auto
	vcftools --gzvcf /data3/hirundo/analysis/ld/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.phased.vcf.gz --keep popmap.rand.$pop --chr NC_053488.1 --maf 0.05 --hap-r2 --ld-window-bp 25000 --out ./r2/r2.$pop.chrZ
done
