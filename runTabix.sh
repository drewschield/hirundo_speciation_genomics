list=$1
for gvcf in `cat $list`; do
	echo indexing $gvcf
	tabix -p vcf $gvcf
done