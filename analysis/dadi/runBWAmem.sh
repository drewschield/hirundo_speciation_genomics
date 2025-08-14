for indv in `cat sample.filter.list`; do
	name=HR${indv}
	echo "Mapping filtered $name data to reference"
	bwa mem -t 16 -R "@RG\tID:$name\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:$name" ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.fasta ./fastq_filtered/${indv}.trim.fq.gz | samtools sort -@ 16 -O bam -T $name.temp -o ./bam/$name.bam -
	samtools index -@ 8 ./bam/$name.bam
done
