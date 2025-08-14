for bam in ./bam/*.bam; do
	name=${bam##*/}	# removes directory name
	echo calculating mapping statistics for $name
	samtools stat -@ 8 $bam > ./analysis/mapping_statistics/$name.stat.txt
done
