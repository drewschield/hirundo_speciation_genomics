for chrom in `cat chrom.list`; do
	Rscript rehhScans.R $chrom
done
