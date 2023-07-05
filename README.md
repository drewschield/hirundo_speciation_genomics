# Sexual selection promotes reproductive isolation in barn swallows: Genomic analysis workflow

![Sexual selection promotes reproductive isolation in barn swallows](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps used to study the role of sexual selection in speciation through 1) the characterization of the genetic architecture of sexual signal traits in barn swallows (_Hirundo rustica_), 2) tests of divergent sexual selection in geographic isolation, 3) tests of selection against gene flow at trait loci in secondary contact, and 4) tests that selection has maintained associations between trait loci, known as genetic coupling, leading to reproductive isolation due to sexual selection. This workflow is a companion to the methods described in Schield et al. (_in review_).

Lists and reference files can be found in various analysis directories, as can Shell and Python scripts used in elements of the analysis workflow. R scripts are in the `R` directory. Note that you may need to adjust the organization of path/file locations to suit your environment. This workflow assumes that you return to the main working directory `~/hirundo_speciation_genomics` after completing major analysis sections. Adjust the path before the main working directory as necessary. Walkthrough of software installations will take place in `~/hirundo_speciation_genomics/tmp/`.

Feel free to contact me at drew.schield[at]colorado.edu with any questions.

## Contents

* [Software and dependencies](#software-and-dependencies)
* [Genome data processing and variant calling](#genome-data-processing-and-variant-calling)
* [Mapping statistics](#mapping-statistics)
* [PCA](#pca)
* [ADMIXTURE](#admixture)
* [Hybrid index and heterozygosity](#hybrid-index-and-heterozygosity)
* [Demographic inference](#demographic-inference)
* [Appendix](#appendix)
	* [Assignment of B10K barn swallow genome scaffolds to chromosomes](#assignment-of-b10k-barn-swallow-genome-scaffolds-to-chromosomes)
	* [Repeat masking the reference genome](#repeat-masking-the-reference-genome)

## Software and dependencies

The analysis sections below use the following software and dependencies and assume they are on the user path unless otherwise specified:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bwa](http://bio-bwa.sourceforge.net/)
* [htslib](http://www.htslib.org/)
* [Samtools](http://www.htslib.org/)
* [bgzip](http://www.htslib.org/)
* [tabix](http://www.htslib.org/)
* [GATK](https://gatk.broadinstitute.org/hc/en-us)
* [Picard](https://broadinstitute.github.io/picard/)
* [bcftools](http://www.htslib.org/)
* [Plink](https://www.cog-genomics.org/plink/)
* [vcftools](https://vcftools.github.io/index.html)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [SNPRelate](http://bioconductor.org/packages/release/bioc/html/SNPRelate.html)
* [ADMIXTURE](https://dalexander.github.io/admixture/publications.html)
* [Introgress](https://www.uwyo.edu/buerkle/software/introgress/)
* [dadi](https://dadi.readthedocs.io/en/latest/)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [rehh](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html)
* [MashMap](https://github.com/marbl/MashMap)
* [Repeatmasker](https://www.repeatmasker.org/)
* [R](https://cran.r-project.org/)

Note, I installed a number of these programs to my [conda](https://docs.conda.io/en/latest/) environment, or installed via a virtual environment (details below).

[Back to top](#contents)

## Genome data processing and variant calling

This section includes details on steps to filter raw fastq data, map reads to the reference genome, call variants among samples, and filter variant calls. Raw read data are available from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA323498), and should be placed in the `fastq` subdirectory below. Additional information about the reference genome and annotation can be found in the [Appendix](#appendix). 

Note: scripts used in this section are in the `scripts` directory (for tidiness) and will need to be moved to the working directory prior to running this part of the workflow.

### Set up environment

We'll set up directories for read data, mapping data, and variant calls. The `log` subdirectory will contain sample lists, log files, etc.
```
mkdir fastq
mkdir fastq_filtered
mkdir bam
mkdir gvcf
mkdir vcf
mkdir log
```

### Filter raw read data with Trimmomatic

We'll use `Trimmomatic` to filter and quality trim raw read data, using the settings:
* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30

We'll process the input data slightly differently for data generated in 2018 and 2021, based on some slight differences in naming conventions between sequencing runs.

#### Input sample lists

Format `./log/trimmomatic.2018.list` and `./log/trimmomatic.2021.list`.

#### Run Trimmomatic

1. Run `runTrimmomatic.2018.sh` to process 2018 fastq data.
```
sh runTrimmomatic.2018.sh > ./log/runTrimmomatic.2018.log
```

2. Run `runTrimmomatic.2021.sh` to process 2021. fastq data.
```
sh runTrimmomatic.2021.sh > runTrimmomatic.2021.log
```

3. Run `trimmomatic` on outgroup H. smithii individual.
```
trimmomatic PE -phred33 -threads 16 ./fastq/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1.fq.gz ./fastq/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1_P.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1_U.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2_P.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
```

#### Remove unpaired filtered reads

We won't map these.
```
cd fastq_filtered
rm ./*U.trim.fq.gz
cd ..
```

Analysis-read paired read data are in `fastq_filtered`.

### Map reads using BWA

We'll run BWA `mem` on each sample to map it to the reference genome.

#### Run BWA on ingroup samples

1. Format `./log/bwa_mem.list` with all ingroup samples.

2. The script `runBWAmem.sh` runs BWA `mem` and also indexes using `samtools`.
```
sh runBWAmem.sh > ./log/runBWAmem.log
```

#### Run BWA on outgroup H. smithii.

```
bwa mem -t 16 -R "@RG\tID:RS_5\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:RS_5" Hirundo_rustica_bHirRus1.final.fasta ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_1_P.trim.fq.gz ./fastq_filtered/RS_5_USPD16098865-AK1579-AK402_HL7CWDSXX_L1_2_P.trim.fq.gz | samtools sort -@ 16 -O bam -T RS_5.temp -o ./bam/RS_5.bam -
samtools index ./bam/RS_5.bam
```

Analysis-ready bam files are in `./bam`.

### Call raw variants using GATK

We'll call individual variants using `HaplotypeCaller` and cohort variants using `GenotypeGVCFs`.

#### Call individual variants

Run runGATKhaplotypecaller.sh to run `HaplotypeCaller` on samples in `./log/GATKhaplotypecaller.list`.
```
sh runGATKhaplotypecaller.sh ./log/GATKhaplotypecaller.list > ./log/runGATKhaplotypecaller.log
```

The output genomic VCFs for each individual are in `./gvcf/`.

Note, this takes a small eternity to run on all samples in sequence. We can split out the job to parallel processes by breaking up the list of samples and running the script on each list separately (constrained by CPU threads).

#### Call cohort variants

1. Format `./log/sample.hrustica+smithii.gvcf.list`, with the list of paths to all input gvcf files.

2. Run `runTabix.sh` to index all input gvcf files.
```
sh tabix.sh ./log/sample.hrustica+smithii.gvcf.list > ./log/tabix.log
```

3. Format intervals files with lists of scaffolds to run parallel analyses on in `./log/scaffold.list1.intervals`, etc. (13 files total). These are split out based on the largest scaffold, which represents ~14% of the genome.

4. Run GATK `GenotypeGVCFs` to call cohort variants on the scaffold intervals.
```
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list1.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list1.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list1.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list2.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list2.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list3.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list3.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list3.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list4.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list4.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list4.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list5.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list5.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list5.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list6.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list6.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list6.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list7.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list7.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list7.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list8.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list8.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list8.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list9.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list9.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list9.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list10.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list10.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list10.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list11.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list11.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list11.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list12.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list12.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list12.vcf.log
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -L ./log/scaffold.list13.intervals -V ./sample.hrustica+smithii.gvcf.list -allSites -o ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list13.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list13.vcf.log
```

5. Check on running progress.
```
for i in ./log/GenotypeGVCFs.hirundo_rustica+smithii.allsites.raw.scaffold.list*.log; do echo $i; grep 'INFO' $i | tail -n 1; done
```

### Variant filtering

We'll filter to retain a set of high-quality SNP calls, while also maintaining an 'all-sites' VCF through most steps, enabling calculation of statistics requiring both variant and invariant genotypes. We'll also apply several downstream filters for specific analyses (e.g., filtering by minor allele frequency or thinning by physical distance, etc.).

#### Hard filters

We have raw variant calls in 13 separate interval files from the step above. We'll also set hard filters on genotypes in parallel to reduce runtime, then merge the VCFs.

1. Run GATK `VariantFiltration` to set hard filters.
```
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list1.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list1.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list1.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list2.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list2.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list2.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list3.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list3.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list3.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list4.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list4.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list4.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list5.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list5.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list5.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list6.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list6.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list6.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list7.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list7.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list7.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list8.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list8.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list8.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list9.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list9.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list9.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list10.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list10.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list10.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list11.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list11.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list11.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list12.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list12.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list12.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list13.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.list13.vcf.gz > ./log/GATK_VariantFiltration.scaffold.list13.log
``` 

2. Format `./log/vcf.HardFilter.list` to direct to files to be merged.

3. Run Picard `MergeVCFs` to merge the VCF files.
```
java -jar picard.jar MergeVcfs -I ./log/vcf.HardFilter.list -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.vcf.gz
```

4. Remove raw interval VCFs taking up disk space.
```
rm ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.list*
```

5. Mask indels and hard filtered genotypes.
```
bcftools filter --threads 20 -e 'TYPE="indel" || FILTER="QD2" || FILTER="FS60" || FILTER="MQ40" || FILTER="MQRankSum-12.5" || FILTER="ReadPosRankSum-8"' --set-GTs . -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.vcf.gz
bcftools index --threads 48 -t ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz
```

#### Filter samples with high amounts of missing data

1. Query a subset of the filtered VCF to quantify proportion of missing genotypes across samples.
```
bcftools view --threads 16 -r NC_053451.1 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.snps.chr1.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz > ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.snps.chr1.log
vcftools --gzvcf ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.snps.chr1.vcf.gz --missing-indv --out ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.snps.chr1
rm ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.snps.chr1.vcf.gz
```

Based on this analysis, samples HRUWBM90116 and HRUWBM87931 had > 75% missing data.

2. Format `./log/sample.remove.list`.

3. Remove samples in list.
```
bcftools view --threads 16 -S ^./log/sample.remove.list -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.vcf.gz
```

#### Processing sex chromosomes

Females are hemizygous for the Z chromosome, and so cannot have heterozygous sites on Z-linked scaffolds. We'll extract SNPs for Z-linked scaffolds, then identify SNPs in which any known females have heterozygous genotypes, then mask these sites in all individuals. We also want to remove males from the W chromosome, as they don't have one, and also to mask heterozygous sites on the W in females, since these can't exist.

We'll begin by extracting VCF tables for chromosome-assigned autosome, Z, and W scaffolds, respectively. From here forward we'll also omit un-assigned scaffolds from analysis.

1. Format autosome, Z chromosome, and W chromosome-assigned scaffold BED files.
```
./log/scaffold.assigned-autosome.bed
./log/scaffold.assigned-chrZ.bed
./log/scaffold.assigned-chrW.bed
```

2. Extract VCFs based on BED files.
```
bcftools view --threads 16 -R ./log/scaffold.assigned-autosome.bed -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.auto.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.vcf.gz
bcftools view --threads 16 -R ./log/scaffold.assigned-chrZ.bed -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.vcf.gz
bcftools view --threads 16 -R ./log/scaffold.assigned-chrW.bed -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.vcf.gz
```

3. Rename and index autosome all-sites VCF.
```
mv ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.auto.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.final.auto.vcf.gz
```

4. Mask female heterozygous sites on the Z chromosome.

* Extract biallelic SNPs on Z-linked scaffolds.
```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.vcf.gz > ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.snps.chrZ.vcf.gz
```

* Format `./log/sample.female.list`.

* Run `ZchrFemaleHeterozygous.py` to identify and write positions of the Z chromosome with heterozygous genotype calls in any females.
```
gunzip ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.vcf.gz
python ZchrFemaleHeterozygous.py ./log/sample.female.list ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.vcf ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.FemaleZhetSites.txt
```

* Convert output to BED format and index it as a GATK feature file.
```
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.FemaleZhetSites.txt > ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.FemaleZhetSites.bed
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.FemaleZhetSites.bed
```

* Run GATK `VariantFiltration` to annotate female heterozygous sites and mask with `bcftools`.
```
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.vcf.gz
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.vcf.gz --mask ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.snps.FemaleZhetSites.bed --mask-name ZHET -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.filter.vcf.gz > ./log/GATK_VariantFiltration_hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.filter.log
bcftools filter --threads 16 -e 'FILTER="ZHET"' --set-GTs . -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrZ.filter.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz
```

5. Mask female heterozygous sites on the W chromosome.

* Remove males from W chromosome VCF.
```
bcftools view --threads 16 -S ./log/sample.female.list -O z -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.vcf.gz
```

* Extract W-linked SNPs.
```
bcftools view --threads 16 -m2 -M2 -U -v snps -O v -o ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.vcf ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.vcf.gz
```

* Identify heterozygous sites on the W chromosome.
```
python ZchrFemaleHeterozygous.py ./log/sample.female.list ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.vcf ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.FemaleHetSites.txt
```

* Convert output to BED format and index it as a GATK feature file.
```
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.FemaleHetSites.txt > ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.FemaleHetSites.bed
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.FemaleHetSites.bed
```

* Run GATK `VariantFiltration` to annotate female heterozygous sites and mask with `bcftools`.
```
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.vcf.gz
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.vcf.gz --mask ./log/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.snps.FemaleHetSites.bed --mask-name WHET -O ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.filter.vcf.gz > ./log/GATK_VariantFiltration_hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.filter.log &
bcftools filter --threads 16 -e 'FILTER="WHET"' --set-GTs . -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.HardFilter.recode.indv.chrW.female.filter.vcf.gz &
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz
```

#### Repeat masking

The steps to annotate repeats using RepeatMaster are detailed in the [Appendix](#appendix).

1. Index repeat BED file.
```
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file ./genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed
```

2. Annotate repeats using GATK `VariantFiltration`.
```
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.final.auto.vcf.gz --mask ./genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O ./vcf/hirundo_rustica+smithii.allsites.final.auto.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_rustica+smithii.auto.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --mask ./genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_rustica+smithii.chrZ.log
./gatk-4.0.8.1/gatk VariantFiltration -V ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz --mask ./genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed --mask-name REP -O ./vcf/hirundo_rustica+smithii.allsites.final.chrW.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.hirundo_rustica+smithii.chrW.log
```

3. Recode repeats as missing genotypes and index VCFs.
```
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto.tmp-rep.vcf.gz
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.tmp-rep.vcf.gz
bcftools filter --threads 24 -e 'FILTER="REP"' --set-GTs . -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.tmp-rep.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.final.auto.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz 
```

#### Post-processing and additional filters for analysis

1. Parse chromosome-specific all-sites VCFs.

* Set up environment
```
cd vcf
mkdir chrom-specific
cd ..
```

* Format chromosome-scaffold table for reference: `./log/chromosome-scaffold.table.txt`

* Run `parseVCFchrW.sh` to parse W chromosome scaffolds.
```
sh parseVCFchrW.sh
```

* Run `parseVCFchrZ.sh` to parse Z chromosome scaffolds.
```
sh parseVCFchrZ.sh
```

* Run `parseVCFauto.sh` to parse autosome scaffolds.
```
sh parseVCFauto.sh
```

2. Extract biallelic SNPs.
```
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto.vcf.gz
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz
bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.vcf.gz
```

3. Extract SNPs with < 20% missing data.
```
bcftools view -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.vcf.gz
bcftools view -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.vcf.gz
bcftools view -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.vcf.gz
```

4. Impose minor allele frequency and minor allele count filters.
```
bcftools view --threads 16 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.vcf.gz
bcftools view --threads 16 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.maf05.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.vcf.gz
bcftools view --threads 16 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.vcf.gz
```

5. Remove H. smithii outgroup.
```
bcftools view --threads 16 -s ^RS_5 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.ingroup.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.vcf.gz
bcftools view --threads 16 -s ^RS_5 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.maf05.ingroup.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.maf05.vcf.gz
bcftools view --threads 16 -s ^RS_5 -c 2 -q 0.05:minor -m2 -M2 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.ingroup.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.vcf.gz
```

#### Cleaning up

These filtering steps produce many intermediate VCF files, especially in the initial parallel processing steps for raw and hard filtered VCFs. Delete these to save some disk space.
```
rm ./vcf/hirundo_rustica+smithii.allsites.raw.scaffold.*
rm ./vcf/hirundo_rustica+smithii.allsites.HardFilter.scaffold.*
```

[Back to top](#contents)

## Mapping statistics

### Set up environment
```
mkdir analysis
cd analysis
mkdir mapping_statistics
cd ..
```

### Calculate mapping statistics

Run `samtoolsStat.sh` from the `scripts` directory. The results are written to `./analysis/mapping_statistics`.
```
sh ./scripts/samtoolsStat.sh 
```

### Parse number of bases mapped per sample

```
for i in ./analysis/mapping_statistics/*.stat.txt; do name=`echo $i | cut -d. -f1`; bases=`grep 'bases mapped (cigar):' $i | cut -f 3`; echo -e "$name\t$bases"; done
```

[Back to top](#contents)

## PCA

We'll examine population genetic structure within barn swallows using the R package `SNPRelate`.

### Set up environment
```
cd ./analysis
mkdir pca
cd pca
```

### Filter input VCF

We'll run PCA on autosomal and Z-linked SNPs after performing several filtering steps:
* Remove high missing data samples
* Remove mislabeled individuals from the sequencing step
* Retail SNPs with a minor allele frequency of >= 0.1

1. Merge autosomal and Z-linked VCFs.
```
bcftools concat -O z -o hirundo_rustica+smithii.allsites.final.snps.miss02.maf05.ingroup.auto+chrZ.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.ingroup.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.maf05.ingroup.vcf.gz
```

2. Remove samples, filter by minor allele frequency.

There were a small number of _savignii_ and _erythrogaster_ samples (6) that were mislabeled during sequencing library preparation, making it unclear which subspecies they belonged to. There was also an inflection point in the proportion of missing genotypes within samples at ~40%. We'll remove the mislabeled and high missing data samples from further analysis.

Format `sample.remove.list`.

Remove samples in the list and filter to retain SNPs with minor allele frequency >= 0.1.
```
bcftools view --threads 16 -S ^sample.remove.list -c 2 -q 0.1:minor -m2 -M2 -O z -o hirundo_rustica+smithii.allsites.final.snps.miss02.maf1.ingroup.indv.auto+chrZ.vcf.gz hirundo_rustica+smithii.allsites.final.snps.miss02.maf05.ingroup.auto+chrZ.vcf.gz
```

3. Sample a million SNPs at random for analysis.
```
bcftools query -f '%CHROM\t%POS\n' hirundo_rustica+smithii.allsites.final.snps.miss02.maf1.ingroup.indv.auto+chrZ.vcf.gz | shuf -n 1000000 > sites.subset.list
```

There were actually only 837,275 SNPs after applying the allele frequency filter. Plenty for analysis.

### Perform PCA

We'll run the PCA in R using the script `./R/pca.R`.

[Back to top](#contents)

## ADMIXTURE

We'll estimate individual admixture proportions between one or more genetic clusters using `admixture`.

### Set up environment
```
cd ./analysis/
mkdir admixture
cd admixture
mkdir input
```

### Convert VCF to Plink input format
```
plink --vcf ../pca/hirundo_rustica+smithii.allsites.final.snps.miss02.maf1.ingroup.indv.auto+chrZ.vcf.gz --make-bed --out hirundo.auto+chrZ --allow-extra-chr --recode12
```

### Fix formatting for ADMIXTURE

The program can't take non-integer scaffold names.
```
sed -i.bak -e 's/NC_0//g' ./input/hirundo.auto+chrZ.bim
sed -i.bak -e 's/NW_0//g' ./input/hirundo.auto+chrZ.bim
sed -i.bak -e 's/\.1//g' ./input/hirundo.auto+chrZ.bim
sed -i.bak -e 's/NC_0//g' ./input/hirundo.auto+chrZ.map
sed -i.bak -e 's/NW_0//g' ./input/hirundo.auto+chrZ.map
sed -i.bak -e 's/\.1//g' ./input/hirundo.auto+chrZ.map
```

### Perform analysis

Run `runAdmixture.sh` to run ADMIXTURE across a series of K values 1-10.
```
sh runAdmixture.sh ./input/hirundo.auto+chrZ.ped
```

### Evaluate K model with lowest cross-validation error.
```
grep -h CV log*.out
```

[Back to top](#contents)

## Hybrid index and heterozygosity

We'll compare individual hybrid index with interspecific heterozygosity at ancestry-informative sites to characterize our current sampling of hybrid zone populations. Elsie Shogren kindly shared her R Markdown detailing her procedure for doing this analysis in her system (thank you, Elsie!).

### Set up environment
```
cd ./analysis
mkdir introgress
cd introgress
mkdir fst
mkdir input
```

### Calculate SNP-based Fst between parental populations

We'll use Fst to determine which ancestry-informative SNPs to include in analysis.

1. Format parental, hybrid, and parental+hybrid population maps.
```
popmap.gutturalis
popmap.rustica
popmap.tytleri
popmap.rustica-tytleri
popmap.rustica-gutturalis
popmap.tytleri-gutturalis
popmap.rustica-hybrids-tytleri
popmap.rustica-hybrids-gutturalis
popmap.tytleri-hybrids-gutturalis
```

2. Merge autosomal and Z chromosome VCFs for biallelic SNPs.
```
bcftools concat --threads 24 -O z -o ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.ingroup.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.maf05.ingroup.vcf.gz
```

3. Calculate Fst between parental populations.
```
vcftools --gzvcf ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz --weir-fst-pop popmap.rustica --weir-fst-pop popmap.tytleri --out ./fst/fst_rustica-tytleri
vcftools --gzvcf ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz --weir-fst-pop popmap.rustica --weir-fst-pop popmap.gutturalis --out ./fst/fst_rustica-gutturalis
vcftools --gzvcf ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz --weir-fst-pop popmap.tytleri --weir-fst-pop popmap.gutturalis --out ./fst/fst_tytleri-gutturalis
```

### Identify ancestry-informative SNPs

Here, we'll look for SNPs with Fst >= 0.6 between parental populations. Barn swallows have very few fixed differences, as shown in previous studies (e.g., Safran et al. 2016; Scordato et al. 2017; Schield et al. 2021), but using high-Fst SNPs should provide sufficient information.

1. Format lists of high-Fst SNPs.
```
tail -n+2 ./fst/fst_rustica-tytleri.weir.fst | awk 'BEGIN{OFS="\t"}{if ($3 >=0.6) print $1,$2}' > ./fst/anc-info.snps.rustica-tytleri.txt
tail -n+2 ./fst/fst_rustica-gutturalis.weir.fst | awk 'BEGIN{OFS="\t"}{if ($3 >=0.6) print $1,$2}' > ./fst/anc-info.snps.rustica-gutturalis.txt
tail -n+2 ./fst/fst_tytleri-gutturalis.weir.fst | awk 'BEGIN{OFS="\t"}{if ($3 >=0.6) print $1,$2}' > ./fst/anc-info.snps.tytleri-gutturalis.txt
```

2. Remove SNPs with completely missing data in hybrids.

Pilot analyses revealed that there are a small number of SNPs with completely missing genotypes in hybrid tytleri-gutturalis, which throws an error in downstream `introgress` analysis.

Identify which SNPs are causing the problem.
```
bcftools view -i 'F_MISSING=0' ./input/anc-info.snps.tytleri-gutturalis_hybrids.vcf.gz | bcftools query -f '%CHROM\t%POS\n'
```

This outputs:
```
NC_053488.1	16368106
NC_053488.1	16368109
NC_053488.1	16368122
```

Note, the queried file in the command above was generated after extracting the SNP set for the tytleri-gutturalis hybrids.

Remove these SNPs from the parental input file.
```
grep -vwE "(16368106|16368109|16368122)" ./fst/anc-info.snps.tytleri-gutturalis.txt > ./fst/tmp.anc-info.snps.tytleri-gutturalis.txt
rm ./fst/anc-info.snps.tytleri-gutturalis.txt
mv ./fst/tmp.anc-info.snps.tytleri-gutturalis.txt ./fst/anc-info.snps.tytleri-gutturalis.txt
```

3. Extract ancestry-informative SNPs for parentals and hybrids.
```
bcftools view --threads 16 -S popmap.rustica -R ./fst/anc-info.snps.rustica-tytleri.txt -O z -o ./input/anc-info.snps.rustica-tytleri_rustica.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.tytleri -R ./fst/anc-info.snps.rustica-tytleri.txt -O z -o ./input/anc-info.snps.rustica-tytleri_tytleri.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.rustica-tytleri -R ./fst/anc-info.snps.rustica-tytleri.txt -O z -o ./input/anc-info.snps.rustica-tytleri_hybrids.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.rustica -R ./fst/anc-info.snps.rustica-gutturalis.txt -O z -o ./input/anc-info.snps.rustica-gutturalis_rustica.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.gutturalis -R ./fst/anc-info.snps.rustica-gutturalis.txt -O z -o ./input/anc-info.snps.rustica-gutturalis_gutturalis.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.rustica-gutturalis -R ./fst/anc-info.snps.rustica-gutturalis.txt -O z -o ./input/anc-info.snps.rustica-gutturalis_hybrids.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.tytleri -R ./fst/anc-info.snps.tytleri-gutturalis.txt -O z -o ./input/anc-info.snps.tytleri-gutturalis_tytleri.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.gutturalis -R ./fst/anc-info.snps.tytleri-gutturalis.txt -O z -o ./input/anc-info.snps.tytleri-gutturalis_gutturalis.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
bcftools view --threads 16 -S popmap.tytleri-gutturalis -R ./fst/anc-info.snps.tytleri-gutturalis.txt -O z -o ./input/anc-info.snps.tytleri-gutturalis_hybrids.vcf.gz /media/drewschield/VernalBucket/hirundo/vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz
```

### Estimate hybrid index and interspecific heterozygosity

We'll run the analysis in R using the script `./R/introgress.R`.

[Back to top](#contents)

## Demographic inference

We'll infer two-population demographic histories and the relative timing of secondary contact between parental populations using modifications to `dadi` described in [Rougemont et al. 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13664) and available via [GitHub](https://github.com/QuentinRougemont/DemographicInference). Preliminary analyses revealed that our sample sizes from the whole genome dataset were insufficient to resolve demographic histories following very recent divergence (consistent with simulation-based results from [Robinson et al. 2014](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-014-0254-4)). To capitalize on a larger sample size, we'll reanalyze RADseq data from [Scordato et al. 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14276) and [Scordato et al. 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13420), which were demutiplexed and quality trimmed as described in [Schield et al. 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15885).

### Set up environment
```
cd ./analysis
mkdir dadi
cd dadi
mkdir fastq_filtered
mkdir bam
mkdir gvcf
mkdir vcf
mkdir log
```

### Process RAD data and call variants with outgroup

The filtered fastq data for the RADseq data from previous studies need to be in `./fastq_filtered`. We already have a gVCF for the H. smithii outgroup, which we will incorporate later.

#### Map RAD data to the reference genome

1. Format `sample.filter.popmap`.

2. Format `sample.filter.list`.
```
cut -f1 sample.filter.popmap > sample.filter.list
```

3. Run `runBWAmem.sh` to map to the reference.
```
sh runBWAmem.sh > ./log/runBWAmem.log
```

4. Format `./bam/bam.list` file.

#### Call and filter variants

1. Temporarily increase the limit on number of open files to accommodate bam files.
```
ulimit -n 1000
```

2. Call variants using bcftools `mpileup` and `call`.
```
cd ./bam
bcftools mpileup --threads 16 -Ou -B -a "DP,AD" -b bam.list -f ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.fasta | bcftools call --threads 16 -c -v -V indels -f "GQ" -Oz -o ../vcf/hirundo_rustica.rad.snps.tmp.vcf.gz
tabix -p vcf ../hirundo_rustica.rad.snps.tmp.vcf.gz
cd ..
```

Note: the '-a "DP,AD"' flag is very important; `mpileup` does not output individual per-base  depth statistics by default, and we will use a DP filter downstream to recode and cull missing data among samples, and ultimately to filter SNPs with high proportions of missing data.

3. Format `scaff.list` with list of scaffolds to retain.

4. Run bcftools `+setGT` to set low depth SNPs as missing.
```
bcftools +setGT ./vcf/hirundo_rustica.rad.snps.tmp.vcf.gz -- -t q -n . -i 'FMT/DP<5' | bcftools view --threads 16 -O z -o ./vcf/hirundo_rustica.rad.snps.tmp.missing.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.rad.snps.tmp.missing.vcf.gz
```

5. Run bcftools `filter` to filter to ordered autosomes and remove SNPs with > 20% missing data.
```
bcftools view --threads 16 -R scaff.list -m2 -M2 -U -v snps -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica.rad.snps.auto.miss02.vcf.gz ./vcf/hirundo_rustica.rad.snps.tmp.missing.vcf.gz
tabix -p vcf ./vcf/hirundo_rustica.rad.snps.auto.miss02.vcf.gz
```

6. Extract SNP data for H. smithii.
```
bcftools view --threads 16 -s RS_5 -O z -o ./vcf/hirundo_smithii.auto.snps.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.vcf.gz
tabix -p vcf ./vcf/hirundo_smithii.auto.snps.vcf.gz 
```

7. Format list of paths to VCFs to merge `vcf.merge.list`.

8. Run bcftools `merge` to combine ingroup and outgroup VCFs and filter to retain SNPs with < 20% missing data.
```
bcftools merge --threads 16 -l vcf.merge.list | bcftools view --threads 16 -m2 -M2 -U -v snps -i 'F_MISSING<0.2' -O z -o ./vcf/hirundo_rustica+smithii.rad.snps.tmp.vcf.gz
```

9. Write SNP positions where H. smithii has missing data to a file.
```
bcftools view -H -s RS_5 ./vcf/hirundo_rustica+smithii.rad.snps.tmp.vcf.gz | grep -v "0/1" | grep -v "0/0" | grep -v "1/1" | cut -d$'\t' -f1,2 > snp.remove.list
```

10. Remove SNPs with missing genotypes for outgroup H. smithii and thin to retain a single SNP per RAD locus.
```
vcftools --gzvcf ./vcf/hirundo_rustica+smithii.rad.snps.tmp.vcf.gz --exclude-positions snp.remove.list --thin 10000 --recode --stdout | bgzip -c > ./vcf/hirundo_rustica+smithii.rad.snps.dadi.vcf.gz
```

11. Write a simplified VCF (to be read into polarization script).
```
(bcftools view -h ./vcf/hirundo_rustica+smithii.rad.snps.dadi.vcf.gz; bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\tGT\\t[%GT\\t]\\n" ./vcf/hirundo_rustica+smithii.rad.snps.dadi.vcf.gz) | cat | bcftools view -m2 -M2 -v snps -O z -o ./vcf/hirundo_rustica+smithii.rad.snps.dadi.fix.vcf.gz
```

12. Polarize VCF by outgroup using the script written by Kristian Ullrich available [here](https://github.com/kullrich/bio-scripts/blob/master/vcf/polarizeVCFbyOutgroup.py).
```
python polarizeVCFbyOutgroup.py -vcf ./vcf/hirundo_rustica+smithii.rad.snps.dadi.fix.vcf.gz -out ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -ind 1 -add
```



































## Appendix

### Assignment of B10K barn swallow genome scaffolds to chromosomes

The [B10K barn swallow reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_015227805.1) includes chromosome-length scaffolds, but which are not assigned to the passerine karyotype. Here, we'll establish synteny between scaffolds and chromosome-assigned scaffolds in the [zebra finch reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_008822105.2/), and reformat the barn swallow scaffolds according to these assignments.

#### Set up environment

```
mkdir reference
cd reference
mkdir reference_B10K
mkdir reference_taeniopygia
cd reference_B10K
mkdir ncbi_version
```

#### Download B10K reference genome

```
cd ncbi_version
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/805/GCF_015227805.1_bHirRus1.pri.v2/GCF_015227805.1_bHirRus1.pri.v2_genomic.fna.gz
gunzip GCF_015227805.1_bHirRus1.pri.v2_genomic.fna.gz
cd ../..
```

#### Download zebra finch genome

```
cd reference_taeniopygia
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna.gz
gunzip GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna.gz
cd ..
```

#### Format zebra finch scaffold data

```
cd reference_taeniopygia
```

Run `fasta_seq_length.py` to calculate scaffold lengths.
```
python fasta_seq_length.py GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna > GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.scaffold_lengths.txt
```

Make list of chromosome-assigned scaffolds `GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrom_scaffolds.list`.

Extract chromosome-assigned scaffolds using `scaffold_list_extractor.py`.
```
python scaffold_list_extractor.py GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrom_scaffolds.list taeniopygia_guttata_ChromAssigned.fasta
cd ..
```

#### Align barn swallow to zebra finch

We'll use MashMap to generate alignments between the reference genomes.

```
cd reference_B10K
mkdir mashmap_results
cd mashmap_results
mashmap -t 4 -r ../../reference_taeniopygia/taeniopygia_guttata_ChromAssigned.fasta -q ../ncbi_version/GCF_015227805.1_bHirRus1.pri.v2_genomic.fna -f one-to-one -s 50000 --pi 90 -o B10K_hirundo_2_zebra_finch.txt
cd ..
```

#### Reformat B10K reference

Format `B10K_scaffold-chromosome_table.txt`.

Format `taeniopygia_scaffold-chromosome_table.txt`

Based on MashMap output, the related scaffolds are in `B10K_barn_swallow_reformat_table.txt`.

Run `reformat_B10K_scaffolds.py` to reformat reference.
```
python reformat_B10K_scaffolds.py B10K_barn_swallow_reformat_table.txt ncbi_version/GCF_015227805.1_bHirRus1.pri.v2_genomic.fna Hirundo_rustica_bHirRus1.final.fasta
cd ..
```

This outputs the original B10K scaffolds, but with correctly-assigned chromosome IDs, and in order of correct chromosome assignments. The chromosome-assigned scaffolds are first, followed by unplaced scaffolds.

#### Set up reference genome for analysis

Here, we'll retrieve the reference to the main working directory and prepare indexes for various analyses (e.g., bwa, GATK, samtools).

```
cd ~/hirundo_speciation_genomics
mv reference/reference_B10K/Hirundo_rustica_bHirRus1.final.fasta .
bwa index Hirundo_rustica_bHirRus1.final.fasta
./gatk-4.0.8.1/gatk CreateSequenceDictionary -R Hirundo_rustica_bHirRus1.final.fasta
samtools faidx Hirundo_rustica_bHirRus1.final.fasta
```

#### Set up reference genome annotation files

```
mkdir genome_annotation
cd genome_annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/805/GCF_015227805.1_bHirRus1.pri.v2/GCF_015227805.1_bHirRus1.pri.v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/805/GCF_015227805.1_bHirRus1.pri.v2/GCF_015227805.1_bHirRus1.pri.v2_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/805/GCF_015227805.1_bHirRus1.pri.v2/GCF_015227805.1_bHirRus1.pri.v2_rna.fna.gz
gunzip *.gz
cd ..
```

[Back to top](#contents)

### Repeat masking the reference genome

We want to avoid analyzing genotypes in repetitive regions of the genome. We'll generate a repeat annotation for filtering using `repeatmasker`.

#### Install repeatmasker dependencies

The dependencies are `rmblast`, `tandem repeats finder`, and `h5py`.

rmblast:
```
cd ./tmp/
wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz
tar -xf rmblast-2.11.0+-x64-linux.tar.gz
sudo cp -r rmblast-2.11.0 /usr/local/
```

Download [tandem repeats finder](https://tandem.bu.edu/trf/trf409.linux64.download.html).
```
sudo cp trf409.linux64 /usr/local/
```

h5py:
```
conda create --name h5py python=3.5
conda deactivate
conda activate h5py
conda install h5py
```

#### Install repeatmasker

Set up and unpack distribution.
```
cd ./tmp/
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz
tar -xf RepeatMasker-4.1.2-p1.tar.gz 
sudo cp -r RepeatMasker /usr/local/
cd /usr/local/RepeatMasker
```

Download expanded repeat library.
```
sudo wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz
sudo gunzip Dfam.h5.gz
sudo mv Dfam.h5 ./Libraries
```

Configure.
```
sudo perl ./configure
```

These paths need to be entered during the configuration:
* `/usr/local/trf409.linux64`
* `/usr/local/rmblast-2.11.0/bin`

#### Repeat mask the reference genome

Run repeatmasker.
```
cd /data3/hirundo/genome_annotation/
mkdir repeatmasker
cd repeatmasker
/usr/local/RepeatMasker/RepeatMasker -pa 32 -xsmall -species aves ../../Hirundo_rustica_bHirRus1.final.fasta -dir .
```

Convert .out file to BED format.
```
tail -n+4 Hirundo_rustica_bHirRus1.final.fasta.out | awk 'BEGIN{OFS="\t"}{print $5,$6-1,$7}' | bedtools sort -i - > Hirundo_rustica_bHirRus1.final.fasta.repeat.bed
```

[Back to top](#contents)





































