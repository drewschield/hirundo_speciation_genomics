list=$1
for i in `cat $list`; do
	./gatk-4.0.8.1/gatk HaplotypeCaller -R Hirundo_rustica_bHirRus1.final.fasta --ERC GVCF -I ./bam/$i.bam -O ./gvcf/$i.raw.snps.indels.g.vcf
	bgzip ./gvcf/$i.raw.snps.indels.g.vcf
done