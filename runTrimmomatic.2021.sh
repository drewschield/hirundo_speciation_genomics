for line in `cat ./log/trimmomatic.2021.list`; do
	name=$line
	echo processing and filtering HR${name}.
	trimmomatic PE -phred33 -threads 16 ./fastq/${name}_*_1.fq.gz ./fastq/${name}_*_2.fq.gz ./fastq_filtered/HR${name}_1_P.trim.fq.gz ./fastq_filtered/HR${name}_1_U.trim.fq.gz ./fastq_filtered/HR${name}_2_P.trim.fq.gz ./fastq_filtered/HR${name}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
done
