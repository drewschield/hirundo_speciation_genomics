for line in `cat ./log/trimmomatic.2018.list`; do
	name=$line
	echo processing and filtering HR${name}.
	cat ./fastq/L_${name}_*_1.fq.gz > ./fastq/HR${name}_1.fq.gz
	cat ./fastq/L_${name}_*_2.fq.gz > ./fastq/HR${name}_2.fq.gz
	trimmomatic PE -phred33 -threads 16 ./fastq/HR${name}_1.fq.gz ./fastq/HR${name}_2.fq.gz ./fastq_filtered/HR${name}_1_P.trim.fq.gz ./fastq_filtered/HR${name}_1_U.trim.fq.gz ./fastq_filtered/HR${name}_2_P.trim.fq.gz ./fastq_filtered/HR${name}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
	rm ./fastq/HR${name}_*.fq.gz
done
