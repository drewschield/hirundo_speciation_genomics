bcftools view --threads 16 -r NC_053488.1 -O z -o ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz
tabix -p vcf ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz
bcftools view --threads 16 -r NW_024403838.1 -O z -o ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ-un1.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz
tabix -p vcf ./vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ-un1.vcf.gz

