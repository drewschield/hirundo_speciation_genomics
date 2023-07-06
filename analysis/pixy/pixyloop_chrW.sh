list=$1
for chrom in `cat $list`; do
	pixy --n_cores 8 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.$chrom.vcf.gz --populations popmap.pixy.female --window_size 200000 --output_folder results --output_prefix pixy_${chrom}_200kb
	pixy --n_cores 8 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.$chrom.vcf.gz --populations popmap.pixy.female --window_size 100000 --output_folder results --output_prefix pixy_${chrom}_100kb
	pixy --n_cores 8 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.$chrom.vcf.gz --populations popmap.pixy.female --window_size 25000 --output_folder results --output_prefix pixy_${chrom}_25kb
	pixy --n_cores 8 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.$chrom.vcf.gz --populations popmap.pixy.female --window_size 10000 --output_folder results --output_prefix pixy_${chrom}_10kb
	pixy --n_cores 8 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.$chrom.vcf.gz --populations popmap.pixy.female --window_size 1000 --output_folder results --output_prefix pixy_${chrom}_1kb
done

