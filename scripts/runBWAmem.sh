for line in `cat ./bwa_mem.list`; do
	name=$line
	echo "Mapping filtered $name data to reference"
	bwa mem -t 16 -R "@RG\tID:$name\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:$name" Hirundo_rustica_bHirRus1.final.fasta ./fastq_filtered/${name}_1_P.trim.fq.gz ./fastq_filtered/${name}_2_P.trim.fq.gz | samtools sort -@ 16 -O bam -T $name.temp -o ./bam/$name.bam -
	samtools index -@ 8 ./bam/$name.bam
done