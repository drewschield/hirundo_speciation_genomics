bcftools view --threads 8 -r NC_053487.1 -O z -o ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz
tabix -p vcf ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz
bcftools view --threads 8 -r NW_024403836.1 -O z -o ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrW-un1.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz
tabix -p vcf ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrW-un1.vcf.gz
bcftools view --threads 8 -r NW_024403837.1 -O z -o ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrW-un2.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz
tabix -p vcf ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrW-un2.vcf.gz

