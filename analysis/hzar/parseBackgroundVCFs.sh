while read region; do
	chrom=`echo "$region" | cut -f 1`;
	start=`echo "$region" | cut -f 2`;
	end=`echo "$region" | cut -f 3`;
	bed=`echo "$chrom:$start-$end"`
	bcftools view --threads 8 -S popmap.rustica-tytleri -r $bed -O z -o ./input_background_introgress/rustica-tytleri.locus.$bed.hybrids.vcf.gz ./tmp.snps.maf1.rustica-tytleri_hybrids.vcf.gz
	bcftools view --threads 8 -S popmap.rustica -r $bed -O z -o ./input_background_introgress/rustica-tytleri.locus.$bed.rustica.vcf.gz ./tmp.snps.maf1.rustica-tytleri_hybrids.vcf.gz
	bcftools view --threads 8 -S popmap.tytleri -r $bed -O z -o ./input_background_introgress/rustica-tytleri.locus.$bed.tytleri.vcf.gz ./tmp.snps.maf1.rustica-tytleri_hybrids.vcf.gz
done < background.rustica-tytleri.snps.maf1.region.txt

while read region; do
	chrom=`echo "$region" | cut -f 1`;
	start=`echo "$region" | cut -f 2`;
	end=`echo "$region" | cut -f 3`;
	bed=`echo "$chrom:$start-$end"`
	bcftools view --threads 8 -S popmap.rustica-gutturalis -r $bed -O z -o ./input_background_introgress/rustica-gutturalis.locus.$bed.hybrids.vcf.gz ./tmp.snps.maf1.rustica-gutturalis_hybrids.vcf.gz
	bcftools view --threads 8 -S popmap.rustica -r $bed -O z -o ./input_background_introgress/rustica-gutturalis.locus.$bed.rustica.vcf.gz ./tmp.snps.maf1.rustica-gutturalis_hybrids.vcf.gz
	bcftools view --threads 8 -S popmap.gutturalis -r $bed -O z -o ./input_background_introgress/rustica-gutturalis.locus.$bed.gutturalis.vcf.gz ./tmp.snps.maf1.rustica-gutturalis_hybrids.vcf.gz
done < background.rustica-gutturalis.snps.maf1.region.txt

while read region; do
	chrom=`echo "$region" | cut -f 1`;
	start=`echo "$region" | cut -f 2`;
	end=`echo "$region" | cut -f 3`;
	bed=`echo "$chrom:$start-$end"`
	bcftools view --threads 8 -S popmap.tytleri-gutturalis -r $bed -O z -o ./input_background_introgress/tytleri-gutturalis.locus.$bed.hybrids.vcf.gz ./tmp.snps.maf1.tytleri-gutturalis_hybrids.vcf.gz
	bcftools view --threads 8 -S popmap.tytleri -r $bed -O z -o ./input_background_introgress/tytleri-gutturalis.locus.$bed.tytleri.vcf.gz ./tmp.snps.maf1.tytleri-gutturalis_hybrids.vcf.gz
	bcftools view --threads 8 -S popmap.gutturalis -r $bed -O z -o ./input_background_introgress/tytleri-gutturalis.locus.$bed.gutturalis.vcf.gz ./tmp.snps.maf1.tytleri-gutturalis_hybrids.vcf.gz
done < background.tytleri-gutturalis.snps.maf1.region.txt
