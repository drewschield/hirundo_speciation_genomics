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
* [Genotype-phenotype associations](#genotype-phenotype-associations)
* [Recombination rate](#recombination-rate)
* [Population genetic diversity and differentiation](#population-genetic-diversity-and-differentiation)
* [Population branch statistics](#population-branch-statistics)
* [Tajima's D](#tajimas-d)
* [Haplotype statistics](#haplotype-statistics)
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
* [easySFS](https://github.com/isaacovercast/easySFS)
* [GEMMA](https://github.com/genetics-statistics/GEMMA)
* [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
* [SMC++](https://github.com/popgenmethods/smcpp)
* [pyrho](https://github.com/popgenmethods/pyrho)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [VCF-kit](https://vcf-kit.readthedocs.io/en/latest/)
* [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/)
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

![PCA](/images/pca.png "pca")

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
mkdir dadi_rustica-tytleri
mkdir dadi_rustica-gutturalis
mkdir dadi_tytleri-gutturalis
```

Distribute copies of the Rougemont et al. `dadi` modifications into each of the `dadi_...` subdirectories.
```
git clone https://github.com/QuentinRougemont/DemographicInference.git
mv DemographicInference dadi_rustica-tytleri

git clone https://github.com/QuentinRougemont/DemographicInference.git
mv DemographicInference dadi_rustica-gutturalis

git clone https://github.com/QuentinRougemont/DemographicInference.git
mv DemographicInference dadi_tytleri-gutturalis
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

12. Polarize VCF by outgroup using the script written by K. Ullrich available [here](https://github.com/kullrich/bio-scripts/blob/master/vcf/polarizeVCFbyOutgroup.py).
```
python polarizeVCFbyOutgroup.py -vcf ./vcf/hirundo_rustica+smithii.rad.snps.dadi.fix.vcf.gz -out ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -ind 1 -add
```

### Format joint site frequency spectrum (jSFS)

We'll convert the SNP data in the polarized VCF to an unfolded jSFS for downstream analysis using `easySFS`

#### Format popmaps

```
popmap.rustica-tytleri
popmap.rustica-gutturalis
popmap.tytleri-gutturalis
```

#### Run projection preview to calculate number of segregating sites with various downsampling schemes
```
conda activate easySFS
./easySFS/easySFS.py -i ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -p popmap.rustica-tytleri --unfolded -a --preview
./easySFS/easySFS.py -i ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -p popmap.rustica-gutturalis --unfolded -a --preview
./easySFS/easySFS.py -i ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -p popmap.tytleri-gutturalis --unfolded -a --preview
```

A sample size of 96 gives a reasonably high number of segregating sites per analysis.

#### Run `easySFS` to output jSFS under projections
```
./easySFS/easySFS.py -i ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -p popmap.rustica-tytleri --unfolded -a -o sfs_rustica-tytleri --prefix sfs.rustica-tytleri --proj 96,96 -f
./easySFS/easySFS.py -i ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -p popmap.rustica-gutturalis --unfolded -a -o sfs_rustica-gutturalis --prefix sfs.rustica-gutturalis --proj 96,96 -f
./easySFS/easySFS.py -i ./vcf/hirundo_rustica+smithii.rad.snps.dadi.polarized.vcf.gz -p popmap.tytleri-gutturalis --unfolded -a -o sfs_tytleri-gutturalis --prefix sfs.tytleri-gutturalis --proj 96,96 -f
conda deactivate
```

#### Copy input jSFS to `03-data` directories
```
cp ./sfs_rustica-tytleri/dadi/rustica-tytleri.sfs ./dadi_rustica-tytleri/03-data/
cp ./sfs_rustica-gutturalis/dadi/rustica-gutturalis.sfs ./dadi_rustica-gutturalis/03-data/
cp ./sfs_tytleri-gutturalis/dadi/tytleri-gutturalis.sfs ./dadi_tytleri-gutturalis/03-data/
```

### Set up environment for dadi analysis

Here, we'll set up a virtual environment with the necessary dependencies for dadi to run based on the workflow and scripts from [Rougemont et al. 2017](https://github.com/QuentinRougemont/DemographicInference).

Dependencies:
* Python 2.7
* scipy 0.13 or older
* numpy
* matplotlib
* pylab
* parallel
* libgfortran3
* mkl

Note, this was not the most straightforward to get set up. Specific details below should help avoid unnecessary hangups.

#### Set up installation directory
```
mkdir dadi-install
```

#### Create new virtual environment for dadi and dependencies
```
cd dadi-install
virtualenv -p /usr/bin/python2.7 dadi-env
source dadi-env/bin/activate
cd ..
```

Run `deactivate` to exit virtual environment.

#### Install scipy 0.13.3 from wheel included with dadi modifications
```
cd dadi_rustica-tytleri/01-scripts
pip install scipy-0.13.3-cp27-cp27mu-linux_x86_64.whl
cd ..
```

#### Install numpy, matplotlib, and parallel
```
pip install numpy
pip install matplotlib==2.0.0
sudo apt-get install parallel
```

#### Install libgfortran3 and mkl; format mkl library file
```
sudo apt-get install libgfortran3
pip install mkl
cd dadi-install/dadi-env/lib/
mv libmkl_rt.so.2 libmkl_rt.so
cd ../../../
```

#### Tell dadi where its libraries are

We need to run this during each new terminal instance prior to running analysis:
```
export LD_LIBRARY_PATH=/data3/hirundo/analysis/dadi/dadi-install/dadi-env/lib/:$LD_LIBRARY_PATH
```

#### Make run scripts executable
```
chmod 777 ./dadi_rustica-tytleri/01-scripts/*.sh
chmod 777 ./dadi_rustica-gutturalis/01-scripts/*.sh
chmod 777 ./dadi_tytleri-gutturalis/01-scripts/*.sh
```

#### Run pilot analysis to make sure things are working
```
cd ./analysis/dadi
source dadi-env/bin/activate
export LD_LIBRARY_PATH=/data3/hirundo/analysis/dadi/dadi-install/dadi-env/lib/:$LD_LIBRARY_PATH
cd dadi_rustica-tytleri
./01-scripts/01-run_model_iteration_v2.sh 1 ./03-data/rustica-tytleri.sfs SI unfolded 80
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs SI unfolded 80
```

Both scripts should run without errors and produce output files in `dadi_rustica-tytleri`.

#### Check that all models will run
```
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs SI unfolded 80
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs AM unfolded 80
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs IM unfolded 80
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs SC unfolded 80
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs AM2m unfolded 80
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs IM2m unfolded 80
sh ./01-scripts/01-run_model_iteration_v2.sh 55 ./03-data/rustica-tytleri.sfs SC2m unfolded 80
```

The 'SC' model fails with "ValueError: need more than 7 values to unpack". If we look at the models detailed in `02-modifs_v2/unfolded/modeledemo_new_models.py` and there are two 'SC' models. The second includes more than 7 parameters (we will comment this out below).

#### Edit models and run scripts

1. Make subdirectory to tinker with scripts.
```
mkdir script_edits
cd script_edits
cp 02-modifs_v2/unfolded/modeledemo_new_models.py .
cp 01-scripts/00.run_dadi_parallel_v2.sh .
```

2. Edit `modeledemo_new_models.py` to comment out second 'SC' model.

3. Edit `00.run_dadi_parallel_v2.sh` to perform 20 iterations.

4. Distribute edited scripts to analysis directories.
```
cp modeledemo_new_models.py ../dadi_rustica-tytleri/02-modifs_v2/unfolded/
cp modeledemo_new_models.py ../dadi_rustica-gutturalis/02-modifs_v2/unfolded/
cp modeledemo_new_models.py ../dadi_tytleri-gutturalis/02-modifs_v2/unfolded/
cp 00.run_dadi_parallel_v2.sh ../dadi_rustica-tytleri/01-scripts/
cp 00.run_dadi_parallel_v2.sh ../dadi_rustica-gutturalis/01-scripts/
cp 00.run_dadi_parallel_v2.sh ../dadi_tytleri-gutturalis/01-scripts/
```

We now have an environment set up to run demographic models.

### Run demographic models

A note on grid size: the `dadi` manual suggests setting a grid size to be somewhat larger than the largest dimension of the SFS analyzed (this point is reiterated [here](https://groups.google.com/g/dadi-user/c/2hSq_Tjicso) by R. Gutenkunst. Larger grid sizes provide higher accuracy at the expense of longer runtimes. The analysis scripts here multiply the input grid size by 2.

#### Run models for rustica-tytleri
```
source ./dadi-install/dadi-env/bin/activate
export LD_LIBRARY_PATH=/data3/hirundo/analysis/dadi/dadi-install/dadi-env/lib/:$LD_LIBRARY_PATH
cd dadi_rustica-tytleri
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs SI unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs AM unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs IM unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs SC unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs AM2m unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs IM2m unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs SC2m unfolded 55
cd ..
```

To check on a running model (e.g., 'SI'):
```
for i in ./10-log/SI_*; do tail $i; done
```

#### Run models for rustica-gutturalis
```
source ./dadi-install/dadi-env/bin/activate
export LD_LIBRARY_PATH=/data3/hirundo/analysis/dadi/dadi-install/dadi-env/lib/:$LD_LIBRARY_PATH
cd dadi_rustica-gutturalis
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs SI unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs AM unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs IM unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs SC unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs AM2m unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs IM2m unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs SC2m unfolded 55
cd ..
```

#### Run models for tytleri-gutturalis
```
source ./dadi-install/dadi-env/bin/activate
export LD_LIBRARY_PATH=/data3/hirundo/analysis/dadi/dadi-install/dadi-env/lib/:$LD_LIBRARY_PATH
cd dadi_tytleri-gutturalis
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs SI unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs AM unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs IM unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs SC unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs AM2m unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs IM2m unfolded 55
./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs SC2m unfolded 55
cd ..
```

### Parse best-fitting models

1. Write `modelLL.sh` to `./analysis/dadi/script_edits` and distribute to analysis directories.
```
cd ./analysis/dadi/script_edits
cp modelLL.sh ../dadi_rustica-tytleri
cp modelLL.sh ../dadi_rustica-gutturalis
cp modelLL.sh ../dadi_tytleri-gutturalis
cd ..
```

2. Run script to parse best-fitting model.
```
cd dadi_rustica-tytleri
sh modelLL.sh
cd ..

cd dadi_rustica-gutturalis
sh modelLL.sh
cd ..

cd dadi_tytleri-gutturalis
sh modelLL.sh
cd ..
```

### Run longer optimization runs for best-fitting models

Here, we'll run 20 iterations of the best-fit models from above, but with larger grid size to increase parameter estimation accuracy.

```
source ./dadi-install/dadi-env/bin/activate
export LD_LIBRARY_PATH=/data3/hirundo/analysis/dadi/dadi-install/dadi-env/lib/:$LD_LIBRARY_PATH
cd dadi_rustica-tytleri
nohup ./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-tytleri.sfs SC2m unfolded 70 &
for i in SC2m_*; do model=`echo $i | cut -d'_' -f1`; iter=`echo $i | cut -d'_' -f2`; ll=`grep 'Optimized log-likelihood' $i/$i.txt | tail -n1`; aic=`grep 'AIC' $i/$i.txt | tail -n1`; theta=`grep 'theta' $i/$i.txt | tail -n1`; echo $model $iter $ll $aic $theta | sort -k3; done
cd ..

cd dadi_rustica-gutturalis
nohup ./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/rustica-gutturalis.sfs SC2m unfolded 70 &
for i in SC2m_*; do model=`echo $i | cut -d'_' -f1`; iter=`echo $i | cut -d'_' -f2`; ll=`grep 'Optimized log-likelihood' $i/$i.txt | tail -n1`; aic=`grep 'AIC' $i/$i.txt | tail -n1`; theta=`grep 'theta' $i/$i.txt | tail -n1`; echo $model $iter $ll $aic $theta | sort -k3; done
cd ..

cd dadi_tytleri-gutturalis
nohup ./01-scripts/00.run_dadi_parallel_v2.sh ./03-data/tytleri-gutturalis.sfs SC unfolded 70 &
for i in SC_*; do model=`echo $i | cut -d'_' -f1`; iter=`echo $i | cut -d'_' -f2`; ll=`grep 'Optimized log-likelihood' $i/$i.txt | tail -n1`; aic=`grep 'AIC' $i/$i.txt | tail -n1`; theta=`grep 'theta' $i/$i.txt | tail -n1`; echo $model $iter $ll $aic $theta | sort -k3; done
cd ..
```

[Back to top](#contents)


## Genotype-phenotype associations

We'll map genetic associations with ventral color and tail streamer length variation using Bayesian sparse linear mixed models (BSLMM) and univariate linear mixed models (LMM) in `GEMMA`.

### Set up environment
```
cd ./analysis/
mkdir gemma
cd gemma
mkdir vcf
mkdir vcf_imputed
mkdir phenotype_lists
mkdir input
```

Format popmaps.
```
popmap.all
popmap.all.male
popmap.all.female
```

### Install GEMMA and Beagle
```
cd ~/hirundo_speciation_genomics/tmp/
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz
chmod +x gemma-0.98.4-linux-static-AMD64
mv gemma-0.98.4-linux-static-AMD64 gemma
sudo cp gemma /usr/local/bin/
cd ../analysis/gemma
wget https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar
```

### Format and impute SNP data

`GEMMA` requires complete genotype information, so we will impute missing genotypes using `beagle`.

#### Concatenate autosomal and Z chromosome VCFs
```
bcftools concat -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.ingroup.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.maf05.ingroup.vcf.gz
```

#### Extract SNP data for males
```
bcftools view --threads 16 -S popmap.all.male -O v ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.vcf.gz
```

#### Extract SNP data for females
```
bcftools view --threads 16 -S popmap.all.female -O v ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.female.vcf.gz
bcftools view --threads 16 -S popmap.all.female -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.female.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.ingroup.vcf.gz
bcftools concat --threads 16 -O z -o ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.female.vcf.gz ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.female.vcf.gz
rm ./vcf/hirundo_rustica+smithii.allsites.final.chrW.snps.miss02.maf05.female.vcf.gz 
rm ./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.female.vcf.gz 
```

#### Perform imputation using Beagle
```
java -Xmx96g -jar beagle.28Jun21.220.jar nthreads=16 gt=./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.vcf.gz out=./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute
java -Xmx96g -jar beagle.28Jun21.220.jar nthreads=16 gt=./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.vcf.gz out=./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute
java -Xmx96g -jar beagle.28Jun21.220.jar nthreads=16 gt=./vcf/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.vcf.gz out=./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute
tabix -p vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz
tabix -p vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.vcf.gz
tabix -p vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.vcf.gz
```

### Format sample lists for various analyses

There are cases where samples have incomplete phenotype matrices. `GEMMA` requires complete data, so we'll downsample our imputed genotypes accordingly.

Sample sizes:
* Hybrid ventral color _n_ = 159
* Hybrid tail streamer length _n_ = 151
* Full ventral color _n_ = 305
* Full tail streamer length _n_ = 300
* Male ventral color _n_ = 157
* Male tail streamer length _n_ = 151
* Female ventral color _n_ = 148
* Female tail streamer length _n_ = 149

Lists of samples in these categories are in `./phenotype_lists/`.

### Extract VCFs and convert to Plink format
```
bcftools view --threads 16 -S ./phenotype_lists/list.hybrid.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.hybrid.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.tail.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.full.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.full.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.tail.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.male.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.male.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.tail.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.female.plum -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.plum.vcf.gz
bcftools view --threads 16 -S ./phenotype_lists/list.female.tail -O v ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.vcf.gz | bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.tail.vcf.gz
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.plum.vcf.gz --make-bed --out ./input/gwas.hybrid.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.hybrid_all.impute.tail.vcf.gz --make-bed --out ./input/gwas.hybrid.tail --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.plum.vcf.gz --make-bed --out ./input/gwas.full.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.ingroup.impute.tail.vcf.gz --make-bed --out ./input/gwas.full.tail --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.plum.vcf.gz --make-bed --out ./input/gwas.male.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ.snps.miss02.maf05.male.impute.tail.vcf.gz --make-bed --out ./input/gwas.male.tail --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.plum.vcf.gz --make-bed --out ./input/gwas.female.plum --allow-extra-chr
plink --vcf ./vcf_imputed/hirundo_rustica+smithii.allsites.final.auto+chrZ+chrW.snps.miss02.maf05.female.impute.tail.vcf.gz --make-bed --out ./input/gwas.female.tail --allow-extra-chr
```

### Format input phenotype data

Here, we'll edit the .fam files generated by `plink` to include the phenotype data.

#### Save copies of unedited .fam files
```
cd input
mkdir original_fam
cp *.fam original_fam
cd ..
```

#### Table of phenotypes in edited .fam files

Note: we've added randomized phenotypes to the .fam files for the full dataset.

| Dataset | File | n1     | n2          | n3    | n4   | n5            |
|---------|------|--------|-------------|-------|------|---------------|
| Hybrid  | plum | breast |             |       |      |               |
| Hybrid  | tail | tail   |             |       |      |               |
| Full    | plum | throat | breast      | belly | vent | breast random |
| Full    | tail | tail   | tail random |       |      |               |
| Male    | plum | throat | breast      | belly | vent |               |
| Male    | tail | tail   |             |       |      |               |
| Female  | plum | throat | breast      | belly | vent |               |
| Female  | tail | tail   |             |       |      |               |

### Generate relatedness matrices
```
gemma -bfile ./input/gwas.hybrid.plum -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.hybrid.plum
gemma -bfile ./input/gwas.hybrid.tail -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.hybrid.tail
gemma -bfile ./input/gwas.full.plum -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.full.plum
gemma -bfile ./input/gwas.full.tail -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.full.tail
gemma -bfile ./input/gwas.male.plum -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.male.plum
gemma -bfile ./input/gwas.male.tail -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.male.tail
gemma -bfile ./input/gwas.female.plum -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.female.plum
gemma -bfile ./input/gwas.female.tail -gk 1 -miss 1 -maf 0 -r2 1 -hwe 0 -o ./gemma.female.tail
```

The output files are written to the auto-generated `./output/` subdirectory.

### BSLMMs

To characterize the overall genetic architecture and proportions of variance explained by all SNPs (PVE) and the proportion of variance explained by sparse effects (alleles of large effect; PGE), we'll run Bayesian sparse linear mixed models in `GEMMA`. We'll run 10 independent chains per analysis with 5 million generations of burn-in and 25 million sampling generations.

#### Run 10 chains per trait

1. Run `runBSLMM.sh` to perform BSLMM analysis.
```
cd output
mkdir log
cd ..
sh runBSLMM.sh ./input/gwas.hybrid.plum ./output/gemma.hybrid.plum.cXX.txt 5000000 25000000 1 hybrid-breast-bright > ./output/log/bslmm.hybrid-breast-bright.log
sh runBSLMM.sh ./input/gwas.hybrid.tail ./output/gemma.hybrid.tail.cXX.txt 5000000 25000000 1 hybrid-tail-streamer > ./output/log/bslmm.hybrid-tail-streamer.log
```

2. Summarize results using R scripts `processBSLMM_full-breast-bright.R` and `processBSLMM_full-tail-streamer.R`, adopted from the [script](https://github.com/edegreef/PUMA-resequencing-data/blob/master/gwas/05-gemma_bslmm_summary.R) from K. Delmore and E. de Greef used in [de Greef et al. 2023](https://www.nature.com/articles/s41598-023-29470-7). They take as input a specified file name convention, then parse and combine results from the multiple runs, outputting summary tables for hyperparameters, as well as specific tables for the parameters, including sets of SNPs passing PIP thresholds and an output of sparse effects for all SNPs (i.e., means of all runs).

```
Rscript processBSLMM_hybrid-breast-bright.R
Rscript processBSLMM_hybrid-tail-streamer.R
```

#### Statistical summary

To summarize the posterior distributions of hyperparameters and per-SNP posterior inclusion probabilities, run `./R/gemma_BSLMM.R`.

### Univariate LMMs

#### Run LMMs for hybrid, full (including randomized phenotypes), male, and female datasets for both traits
```
gemma -bfile ./input/gwas.hybrid.plum -k ./output/gemma.hybrid.plum.cXX.txt -lmm 4 -n 1 -o gwas_full_lmm.hybrid-breast-bright
gemma -bfile ./input/gwas.hybrid.tail -k ./output/gemma.hybrid.tail.cXX.txt -lmm 4 -n 1 -o gwas_full_lmm.hybrid-tail-streamer
gemma -bfile ./input/gwas.full.tail -k ./output/gemma.full.tail.cXX.txt -lmm 4 -n 1 -o gwas_full_lmm.full-tail-streamer
gemma -bfile ./input/gwas.full.plum -k ./output/gemma.full.plum.cXX.txt -lmm 4 -n 2 -o gwas_full_lmm.full-breast-bright
gemma -bfile ./input/gwas.full.tail -k ./output/gemma.full.tail.cXX.txt -lmm 4 -n 2 -o gwas_full_lmm.full-tail-streamer-random
gemma -bfile ./input/gwas.full.plum -k ./output/gemma.full.plum.cXX.txt -lmm 4 -n 5 -o gwas_full_lmm.full-breast-bright-random
gemma -bfile ./input/gwas.male.tail -k ./output/gemma.male.tail.cXX.txt -lmm 4 -n 1 -o gwas_full_lmm.male-tail-streamer
gemma -bfile ./input/gwas.male.plum -k ./output/gemma.male.plum.cXX.txt -lmm 4 -n 2 -o gwas_full_lmm.male-breast-bright
gemma -bfile ./input/gwas.female.tail -k ./output/gemma.female.tail.cXX.txt -lmm 4 -n 1 -o gwas_full_lmm.female-tail-streamer
gemma -bfile ./input/gwas.female.plum -k ./output/gemma.female.plum.cXX.txt -lmm 4 -n 2 -o gwas_full_lmm.female-breast-bright
```

#### Prune the results files

The output files for LMMs are large, so we can prune these down considerably by filtering out SNPs with a -log10(Wald P-value) < 2 (primarily for plotting purposes).

```
mkdir lmm_full_prune
cd output
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.hybrid-breast-bright.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.hybrid-breast-bright.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.hybrid-tail-streamer.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.hybrid-tail-streamer.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.full-tail-streamer.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.full-tail-streamer.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.full-breast-bright.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.full-breast-bright.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.full-tail-streamer-random.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.full-tail-streamer-random.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.full-breast-bright-random.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.full-breast-bright-random.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.male-tail-streamer.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.male-tail-streamer.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.male-breast-bright.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.male-breast-bright.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.female-tail-streamer.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.female-tail-streamer.assoc.prune.txt
awk '{if (-log($13)/log(10)>=2) print $0}' ./gwas_full_lmm.female-breast-bright.assoc.txt | grep -v 'NW_' > ../lmm_full_prune/gwas_full_lmm.female-breast-bright.assoc.prune.txt
cd ..
```

#### Parse significant SNPs

To enable searching of genome features associated with significant SNPs, we'll parse SNPs with a Wald P-value < 0.05 after Bonferroni correction.

Table of SNPs in each model:
| Trait         | Dataset | SNPs    |
|---------------|---------|---------|
| ventral color | hybrid  | 9033285 |
| ventral color | full    | 9311628 |
| ventral color | male    | 8851265 |
| ventral color | female  | 8757718 |
| tail streamer | hybrid  | 8870504 |
| tail streamer | full    | 9246603 |
| tail streamer | male    | 8743581 |
| tail streamer | female  | 8773176 |

```
mkdir significant
cd output
awk '{if ($13 < (0.05/9033285)) print $0}' ./gwas_full_lmm.hybrid-breast-bright.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.hybrid-breast-bright.assoc.sig.txt
awk '{if ($13 < (0.05/8870504)) print $0}' ./gwas_full_lmm.hybrid-tail-streamer.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.hybrid-tail-streamer.assoc.sig.txt
awk '{if ($13 < (0.05/9311628)) print $0}' ./gwas_full_lmm.full-tail-streamer.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.full-tail-streamer.assoc.sig.txt
awk '{if ($13 < (0.05/9246603)) print $0}' ./gwas_full_lmm.full-breast-bright.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.full-breast-bright.assoc.sig.txt
awk '{if ($13 < (0.05/9311628)) print $0}' ./gwas_full_lmm.full-tail-streamer-random.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.full-tail-streamer-random.assoc.sig.txt
awk '{if ($13 < (0.05/9246603)) print $0}' ./gwas_full_lmm.full-breast-bright-random.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.full-breast-bright-random.assoc.sig.txt
awk '{if ($13 < (0.05/8851265)) print $0}' ./gwas_full_lmm.male-tail-streamer.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.male-tail-streamer.assoc.sig.txt
awk '{if ($13 < (0.05/8743581)) print $0}' ./gwas_full_lmm.male-breast-bright.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.male-breast-bright.assoc.sig.txt
awk '{if ($13 < (0.05/8757718)) print $0}' ./gwas_full_lmm.female-tail-streamer.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.female-tail-streamer.assoc.sig.txt
awk '{if ($13 < (0.05/8773176)) print $0}' ./gwas_full_lmm.female-breast-bright.assoc.txt | grep -v 'NW_' > ../significant/gwas_full_lmm.female-breast-bright.assoc.sig.txt
cd ..
```

#### Query genome annotation with significant GWA SNPs
```
mkdir annotation
mkdir significant_annotation
cd annotation
```

1. Make 'gene' and 'exon' GFF files from the full genomic GFF annotation in `~/hirundo_speciation_genomics/genome_annotation`.

```
$awk '{if ($3 == "gene") print $0}' ~/hirundo_speciation_genomics/genome_annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gff > ~/hirundo_speciation_genomics/genome_annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gene.gff
$awk '{if ($3 == "exon") print $0}' ~/hirundo_speciation_genomics/genome_annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gff | bedtools sort -i - > ~/hirundo_speciation_genomics/genome_annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.exon.gff
```

2. Convert significant SNP tables to BED format.
```
cd ./significant
for i in *.sig.txt; do iname=${i%.txt}; tail -n+2 $i | awk 'BEGIN{OFS="\t"}{print $1,$3-1,$3,$13}' | bedtools sort -i - > $iname.bed; done
for i in *-random.assoc.sig.txt; do iname=${i%.txt}; tail -n+2 $i | awk 'BEGIN{OFS="\t"}{print $1,$3-1,$3,$13}' | bedtools sort -i - > $iname.bed; done
```

3. Run bedtools `intersect` to find overlaps between significant SNPs and genes.
```
for i in *.sig.bed; do iname=${i%.bed}; echo -e "chrom-snp\tstart-snp\tend-snp\tWaldP\tchrom-gene\tstart-gene\tend-gene\tstrand\tfeature" > ../significant_annotation/$iname.gene-intersect.txt; bedtools intersect -wo -a $i -b ../annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gene.sort.gff | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$8,$9,$11,$13}' >> ../significant_annotation/$iname.gene-intersect.txt; done
for i in *-random.assoc.sig.bed; do iname=${i%.bed}; echo -e "chrom-snp\tstart-snp\tend-snp\tWaldP\tchrom-gene\tstart-gene\tend-gene\tstrand\tfeature" > ../significant_annotation/$iname.gene-intersect.txt; bedtools intersect -wo -a $i -b ../annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gene.sort.gff | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$8,$9,$11,$13}' >> ../significant_annotation/$iname.gene-intersect.txt; done
```

4. Run bedtools `intersect` to find genes within 50 kb of significant SNPs.
```
for i in *.sig.bed; do iname=${i%.bed}; echo -e "chrom-snp\tstart-snp\tend-snp\tWaldP\tchrom-gene\tstart-gene\tend-gene\tstrand\tfeature" > ../significant_annotation/$iname.gene-intersect+50kb.txt; awk '{OFS="\t"}{print $1,$2-50000,$3+50000,$4}' $i | awk '{OFS="\t"}{print($1,$2<0?0:$2,$3,$4)}' | bedtools intersect -wo -a - -b ../annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gene.sort.gff | awk 'BEGIN{OFS="\t"}{print $1,$2+50000,$3-50000,$4,$5,$8,$9,$11,$13}' >> ../significant_annotation/$iname.gene-intersect+50kb.txt; done
for i in *-random.assoc.sig.bed; do iname=${i%.bed}; echo -e "chrom-snp\tstart-snp\tend-snp\tWaldP\tchrom-gene\tstart-gene\tend-gene\tstrand\tfeature" > ../significant_annotation/$iname.gene-intersect+50kb.txt; awk '{OFS="\t"}{print $1,$2-50000,$3+50000,$4}' $i | awk '{OFS="\t"}{print($1,$2<0?0:$2,$3,$4)}' | bedtools intersect -wo -a - -b ../annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.gene.sort.gff | awk 'BEGIN{OFS="\t"}{print $1,$2+50000,$3-50000,$4,$5,$8,$9,$11,$13}' >> ../significant_annotation/$iname.gene-intersect+50kb.txt; done
```

5. Annotation of genic, coding versus noncoding significant SNPs.
```
tail -n+2 gwas_full_lmm.hybrid-breast-bright.assoc.sig.gene-intersect.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | bedtools closest -d -t first -a - -b ~/hirundo_speciation_genomics/genome_annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.exon.gff > gwas_full_lmm.hybrid-breast-bright.assoc.sig.gene-intersect.exon-distance.txt
tail -n+2 gwas_full_lmm.hybrid-tail-streamer.assoc.sig.gene-intersect.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | bedtools closest -d -t first -a - -b ~/hirundo_speciation_genomics/genome_annotation/GCF_015227805.1_bHirRus1.pri.v2_genomic.exon.gff > gwas_full_lmm.hybrid-tail-streamer.assoc.sig.gene-intersect.exon-distance.txt
```

#### Extract results for specific chromosomes
```
mkdir lmm_full_chrom
for gwas in ./output/gwas_full_lmm*; do iname=${gwas%.txt}; fix=`echo $iname | cut -d'/' -f3`; for chrom in chr1A-NC_053453.1 chr2-NC_053450.1 chrZ-NC_053488.1; do scaff=`echo $chrom | cut -d'-' -f2`; head -n1 $gwas > ./lmm_full_chrom/$fix.$chrom.txt; grep -w $scaff $gwas >> ./lmm_full_chrom/$fix.$chrom.txt; done; done
rm ./lmm_full_chrom/*.log.*
```

#### Statistical summary

To summarize and plot the results, run `./R/gemma_LMM.R`, which also uses the `chr_rename.txt` file.

[Back to top](#contents)


## Recombination rate

We will estimate recombination rate variation across the genome using `pyrho`, which makes use of a population size history inferred using `SMC++`.

### Set up environment
```
cd ./analysis/
mkdir smc++
mkdir pyrho
```

### Population size inference in SMC++

We will use information about population history (i.e., population size at epoch times) to inform recombination rate analysis. SMC++ can infer population history from multiple samples and unphased genotypes. It is also capable of inferring population split times when two populations are analyzed together. Here, we will install SMC++, run it on test data to ensure the build is working properly on the system, then perform analysis on the barn swallows, with the goal of providing a population history to downstream recombination rate inference.

#### Install SMC++

`SMC++` relies on a few [specific dependencies](https://github.com/popgenmethods/smcpp#installation-instructions). We'll install everything together in a virtual environment.

1. Set up install directory
```
cd ~/hirundo_speciation_genomics/tmp/
mkdir smc++-install
cd smc++-install
```
2. Install Python 3.8 and libraries
```
sudo apt install python3.8
sudo apt-get install python3.8-dev
```
3. Create and activate virtual environment for `SMC++`.
```
virtualenv -p /usr/bin/python3.8 smc++-env
source smc++-env/bin/activate
```
4. Install library requirements.
```
sudo apt-get install -y python3-dev libgmp-dev libmpfr-dev libgsl0-dev
```
5. Install `SMC++`.
```
pip install git+https://github.com/popgenmethods/smcpp
```
6. Test that the build worked in the virtual environment.
```
smc++ vcf2smc -h
```

#### Test SMC++ on example data
1. Clone repository with example data.
```
cd ./analysis/smc++
mkdir test_example
cd test_example
mkdir out
mkdir analysis
git clone https://github.com/popgenmethods/smcpp
```
2. Convert VCF to SMC format.
```
source ~/hirundo_speciation_genomics/tmp/smc++-env/bin/activate
smc++ vcf2smc ../smcpp/example/example.vcf.gz out/example.smc.gz 1 Pop1:msp_0,msp_1
```
3. Fit the model using `estimate`.
```
smc++ estimate -o analysis/ 1.25e-8 out/example.smc.gz
```
4. Visualize the results using `plot`.
```
smc++ plot plot.pdf analysis/model.final.json -c
```

#### SMC++ analysis on the barn swallow data

We need a representative population for downstream recombination rate estimation. To avoid issues related to population substructure, we'll analyze rustica from Karasuk, Russia (n = 10). We'll perform analysis using information from all autosomes.

```
cd ./analysis/smc++
mkdir out
mkdir analysis
mkdir vcf
```

1. Format `popmap.rustica.karasuk`.
2. Extract random 5 'distinguished' samples to iterate SMC++ conversion over.
```
shuf -n 5 popmap.rustica.karasuk > distinguished.list
```
3. Format `chromosome-scaffold.auto.table.txt` with conversion between chromosome names and scaffold IDs for autosomes.
4. Generate input VCF.
```
bcftools view --threads 16 -S popmap.rustica.karasuk -c 2:minor -m2 -M2 -U -v snps -O z -o ./vcf/rustica.allsites.final.auto.snps.miss02.mac2.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.vcf.gz
tabix -p vcf ./vcf/rustica.allsites.final.auto.snps.miss02.mac2.vcf.gz
```
5. Extract VCF for each autosome.
```
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; bcftools view --threads 16 -r $scaff -S popmap.rustica.karasuk -O z -o ./vcf/rustica.karasuk.allsites.final.auto.snps.miss02.mac2.$chrom.vcf.gz ./vcf/rustica.allsites.final.auto.snps.miss02.mac2.vcf.gz; done < ./chromosome-scaffold.auto.table.txt
for vcf in ./vcf/*.vcf.gz; do tabix -C -p vcf $vcf; done
```
6. Convert VCFs to SMC input format.
```
source ~/hirundo_speciation_genomics/tmp/smc++-env/bin/activate
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; for indv in `cat distinguished.list`; do smc++ vcf2smc -c 50000 ./vcf/rustica.karasuk.allsites.final.auto.snps.miss02.mac2.$chrom.vcf.gz out/rustica.$chrom.$indv.smc.gz $scaff Pop1:HRVN96101,HRVN96107,HRVN96108,HRVN96103,HRVN96104,HRVN96105,HRVN96106,HRVN96300,HRVN96102,HRVN96298 -d $indv $indv; done; done < ./chromosome-scaffold.auto.table.txt
```
7. Fit the model using `estimate`.
```
smc++ estimate --timepoints 1000 200000 -o analysis/ 2.3e-9 out/rustica.*.smc.gz
```
8. Plot/output results table.
```
smc++ plot rustica-SMC.pdf analysis/model.final.json -g 1 -c
```
This writes `rustica-SMC.csv` and `rustica-SMC.pdf`. The .csv file can be used in downstream `pyrho` analysis.

### Recombination rate inference in pyrho

Here, we will install `pyrho`, run it on test data to ensure the build is working properly on the system, then perform analysis on the barn swallows.

#### Install pyrho

`pyrho` relies on several [specific dependencies](https://github.com/popgenmethods/pyrho). We'll install everything within a virtual environment.

1. Set up install directory.
```
cd ~/hirundo_speciation_genomics/tmp/
mkdir pyrho-install
cd pyrho-install
```
2. Create and activate new virtual environment for `pyrho`.
```
virtualenv -p /usr/bin/python3 pyrho-env
source pyrho-env/bin/activate
```
3. Install `ldpop` dependency.
```
git clone https://github.com/popgenmethods/ldpop.git ldpop
pip install ldpop/
```
4. Install `cython`.
```
pip install cython
```
5. Install msprime and libraries.
```
sudo apt-get install python-dev libgsl0-dev
python3 -m pip install msprime --no-binary msprime
```
6. Install `pyrho`.
```
git clone https://github.com/popgenmethods/pyrho.git pyrho
pip install pyrho/
```
7. Check that install was successful.
```
python -m pytest pyrho/tests/tests.py
```

#### Test pyrho on example data

```
cd ./analysis/pyrho/
mkdir test_example
cd test_example
mkdir out
mkdir analysis
git clone https://github.com/popgenmethods/pyrho.git pyrho
```
1. Precompute a lookup table.
```
pyrho make_table -n 20 -N 25 --mu 1.25e-8 --logfile . --outfile ACB_n_20_N_40_lookuptable.hdf --approx --smcpp_file ACB_pop_sizes.csv --decimate_rel_tol 0.1
```
2. Run `hyperparam` to find hyperparameters that are a good fit to input demography.
```
pyrho hyperparam -n 20 --mu 1.25e-8 --blockpenalty 50,100 --windowsize 25,50 --logfile . --tablefile ACB_n_20_N_40_lookuptable.hdf --num_sims 3 --smcpp_file ACB_pop_sizes.csv --outfile ACB_hyperparam_results.txt 
```
3. Run `optimize` to estimate fine-scale recombination map.
```
pyrho optimize --tablefile ACB_n_20_N_40_lookuptable.hdf --vcffile ACB_chr_1_subset.vcf.gz --outfile ACB_chr_1_subset.rmap --blockpenalty 50 --windowsize 50 --logfile .
```

#### Pyrho analysis on barn swallow data

Now we're set up to run pyrho on each of the ordered chromosomes in the genome. We'll extract chromosome-specific VCFs for rustica to analyze.

```
cd ./analysis/pyrho
mkdir out
mkdir analysis
mkdir vcf
```
1. Format `chromosome-scaffold.auto.table.txt`.
2. Format `chromosome.list` (including Z chromosome).
3. Extract chromosome-specific VCFs.
```
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; bcftools view --threads 16 -r $scaff -O z -o ./vcf/rustica.allsites.final.$chrom.snps.miss02.mac2.vcf.gz ../smc++/vcf/rustica.allsites.final.auto.snps.miss02.mac2.vcf.gz; done < ./chromosome-scaffold.auto.table.txt
bcftools view --threads 16 -S popmap.rustica.karasuk -c 2 -O z -o ./vcf/rustica.allsites.final.chrZ.snps.miss02.mac2.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.vcf.gz
```
4. Generate lookup tables.
```
source ~/hirundo_speciation_genomics/tmp/pyrho-install/pyrho-env/bin/activate
pyrho make_table --numthreads 24 -n 62 -N 78 --mu 2.3e-9 --logfile . --outfile ./lookup/rustica_n_62_N_78_lookuptable.hdf --approx --smcpp_file ../smc++/rustica-SMC.csv
```
5. Find hyperparameter settings that fit the demography.
```
pyrho hyperparam --numthreads 24 -n 62 --mu 2.3e-9 --smcpp_file ../smc++/rustica-SMC.csv --blockpenalty 20,25,50,100 --windowsize 25,50 --logfile . --tablefile ./lookup/rustica_n_62_N_78_lookuptable.hdf --num_sims 3 --outfile ./hyperparam/rustica_hyperparam_results.txt 
```
A window size of 50 and block penalty of 25 look reasonable.
6. Run `optimize` to estimate fine-scale recombination rates.
```
for chrom in `cat chromosome.list`; do pyrho optimize --numthreads 24 --tablefile ./lookup/rustica_n_62_N_78_lookuptable.hdf --vcffile ./vcf/rustica.allsites.final.$chrom.snps.miss02.mac2.vcf.gz --outfile ./results/rustica.$chrom.rmap --blockpenalty 25 --windowsize 50 --ploidy 2 --logfile .; done
```

#### Format results

1. Format window BED files.
```
grep 'NC' ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.auto.genome | bedtools makewindows -g - -w 1000000 -s 100000 > ./Hirundo_rustica_bHirRus1.final.1Mb-100kb.bed; grep 'NC' ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.chrZ.genome | bedtools makewindows -g - -w 1000000 -s 100000 >> ./Hirundo_rustica_bHirRus1.final.1Mb-100kb.bed
grep 'NC' ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.auto.genome | bedtools makewindows -g - -w 1000000 > ./Hirundo_rustica_bHirRus1.final.1Mb.bed; grep 'NC' ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.chrZ.genome | bedtools makewindows -g - -w 1000000 >> ./Hirundo_rustica_bHirRus1.final.1Mb.bed
grep 'NC' ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.auto.genome | bedtools makewindows -g - -w 100000 > ./Hirundo_rustica_bHirRus1.final.100kb.bed; grep 'NC' ~/hirundo_speciation_genomics/Hirundo_rustica_bHirRus1.final.chrZ.genome | bedtools makewindows -g - -w 100000 >> ./Hirundo_rustica_bHirRus1.final.100kb.bed
```
2. Format `chromosome-scaffold.table.txt` (including Z chromosome).
3. Concatenate the recombination map.
```
while read i; do scaff=`echo "$i" | cut -f 1`; chrom=`echo "$i" | cut -f 2`; awk -v var=$scaff 'BEGIN{OFS="\t"}{print var,$1,$2,$3}' ./results/rustica.$chrom.rmap >> ./rustica.all.rmap; done < ./chromosome-scaffold.table.txt
```
4. Calculate mean recombination rate in windows.
```
echo -e "chrom\tstart\tend\trate" > rustica.rmap.1Mb-100kb.txt; bedtools map -a ./Hirundo_rustica_bHirRus1.final.1Mb-100kb.bed -b ./rustica.all.rmap -o mean -c 4 >> rustica.rmap.1Mb-100kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.1Mb.txt; bedtools map -a ./Hirundo_rustica_bHirRus1.final.1Mb.bed -b ./rustica.all.rmap -o mean -c 4 >> rustica.rmap.1Mb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.100kb.txt; bedtools map -a ./Hirundo_rustica_bHirRus1.final.100kb.bed -b ./rustica.all.rmap -o mean -c 4 >> rustica.rmap.100kb.txt
```
5. Format sliding window BED files for specific chromosomes.
```
echo -e 'NC_053453.1\t76187387' | bedtools makewindows -g - -w 100000 -s 10000 > window.100kb-10kb.chr1A-NC_053453.1.bed
echo -e 'NC_053450.1\t156035725' | bedtools makewindows -g - -w 100000 -s 10000 > window.100kb-10kb.chr2-NC_053450.1.bed
echo -e 'NC_053488.1\t90132487' | bedtools makewindows -g - -w 100000 -s 10000 > window.100kb-10kb.chrZ-NC_053488.1.bed
echo -e 'NC_053453.1\t76187387' | bedtools makewindows -g - -w 50000 -s 5000 > window.50kb-5kb.chr1A-NC_053453.1.bed
echo -e 'NC_053450.1\t156035725' | bedtools makewindows -g - -w 50000 -s 5000 > window.50kb-5kb.chr2-NC_053450.1.bed
echo -e 'NC_053488.1\t90132487' | bedtools makewindows -g - -w 50000 -s 5000 > window.50kb-5kb.chrZ-NC_053488.1.bed
echo -e 'NC_053453.1\t76187387' | bedtools makewindows -g - -w 10000 -s 1000 > window.10kb-1kb.chr1A-NC_053453.1.bed
echo -e 'NC_053450.1\t156035725' | bedtools makewindows -g - -w 10000 -s 1000 > window.10kb-1kb.chr2-NC_053450.1.bed
echo -e 'NC_053488.1\t90132487' | bedtools makewindows -g - -w 10000 -s 1000 > window.10kb-1kb.chrZ-NC_053488.1.bed
```
6. Calculate mean recombination rate in sliding windows on specific chromosomes.
```
bedtools sort -i rustica.all.rmap > rustica.all.sort.rmap
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chr1A-NC_053453.1.100kb-10kb.txt; bedtools map -a window.100kb-10kb.chr1A-NC_053453.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chr1A-NC_053453.1.100kb-10kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chr2-NC_053450.1.100kb-10kb.txt; bedtools map -a window.100kb-10kb.chr2-NC_053450.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chr2-NC_053450.1.100kb-10kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chrZ-NC_053488.1.100kb-10kb.txt; bedtools map -a window.100kb-10kb.chrZ-NC_053488.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chrZ-NC_053488.1.100kb-10kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chr1A-NC_053453.1.50kb-5kb.txt; bedtools map -a window.50kb-5kb.chr1A-NC_053453.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chr1A-NC_053453.1.50kb-5kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chr2-NC_053450.1.50kb-5kb.txt; bedtools map -a window.50kb-5kb.chr2-NC_053450.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chr2-NC_053450.1.50kb-5kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chrZ-NC_053488.1.50kb-5kb.txt; bedtools map -a window.50kb-5kb.chrZ-NC_053488.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chrZ-NC_053488.1.50kb-5kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chr1A-NC_053453.1.10kb-1kb.txt; bedtools map -a window.10kb-1kb.chr1A-NC_053453.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chr1A-NC_053453.1.10kb-1kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chr2-NC_053450.1.10kb-1kb.txt; bedtools map -a window.10kb-1kb.chr2-NC_053450.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chr2-NC_053450.1.10kb-1kb.txt
echo -e "chrom\tstart\tend\trate" > rustica.rmap.chrZ-NC_053488.1.10kb-1kb.txt; bedtools map -a window.10kb-1kb.chrZ-NC_053488.1.bed -b rustica.all.sort.rmap -o mean -c 4 >> rustica.rmap.chrZ-NC_053488.1.10kb-1kb.txt
```

#### Summarize results

Run `./R/pyrho.R` to summarize and plot variation in genome-wide recombination rate.

[Back to top](#contents)


## Population genetic diversity and differentiation

We'll use `pixy` to calculate π, Fst, and dxy within and between populations.

### Set up environment
```
cd ./analysis/
mkdir pixy
cd pixy
mkdir results
mkdir results_chrom-specific
```

### Install new version of pixy

A new version of `pixy` was released, with improved run-times and handling of variant data (including bypassing large intermediate zarr files during processing). Runs can also be parallelized. `htslib` is now a dependency.

#### Activate pixy conda environment

We have previously set up `pixy` in its own conda environment with Python 3.6.
```
conda activate pixy
```

#### Remove older version of pixy
```
conda remove pixy
```

#### Install new version
```
conda install -c conda-forge pixy
```

#### Install htslib in environment
```
conda install -c bioconda htslib
```

### Set up population maps

`pixy` uses a population map as a companion input file for summary statistic calculations, which is a two-column tab-delimited file with sample ID and population ID.

Population abbreviations are:

RU = rustica
GU = gutturalis
TY = tytleri
RT = rustica-tytleri
RG = rustica-gutturalis
TG = tytleri-gutturalis
SA = savignii
TR = transitiva
ER = erythrogaster

Population maps are in `popmap.pixy` and `popmap.pixy.subspecies`. Maps for only female samples (for analysis of the W chromosome) are in `popmap.pixy.female` and `popmap.pixy.subspecies.female`.

### Perform analysis

We'll run `pixy` with the chromosome-specific all-sites VCFs in `~/hirundo_speciation_genomics/vcf/chrom-specific` as input.

#### Analysis on autosomes and the Z chromosome

1. Format `scaffold.auto+chrZ.list`.
2. Run `pixyloop.sh` to calculate statistics in windows of various lengths.
```
conda activate pixy
sh pixyloop.sh scaffold.auto+chrZ.list
```
3. Run `pixyloop_subspecies.sh` to calculate statistics in windows with the expanded 'subspecies' popmap.
```
sh pixyloop_subspecies.sh scaffold.auto+chrZ.list
```
#### Analysis on the W chromosome

1. Format `scaffold.chrW.list`.
2. Run `pixyloop_chrW.sh`.
```
sh pixyloop_chrW.sh scaffold.chrW.list
```
2. Run `pixyloop_chrW_subspecies.sh`.
```
sh pixyloop_chrW_subspecies.sh scaffold.chrW.list
```

#### Sliding windows on specific chromosomes

Here we'll run analyses with sliding windows and intermediate step sizes on chromosomes 1A, 2, and the Z chromosome.

|Chromosome |Scaffold     |Length      |
|-----------|-------------|------------|
| 1A        | NC_053453.1 | 76187387   |
| 2         | NC_053450.1 | 156035725  |
| Z         | NC_053488.1 | 90132487   |

1. Format sliding window files for analysis using bedtools `makewindows`.
```
echo -e 'NC_053453.1\t76187387' | bedtools makewindows -g - -w 100000 -s 10000 > window.100kb-10kb.chr1A-NC_053453.1.bed
echo -e 'NC_053450.1\t156035725' | bedtools makewindows -g - -w 100000 -s 10000 > window.100kb-10kb.chr2-NC_053450.1.bed
echo -e 'NC_053488.1\t90132487' | bedtools makewindows -g - -w 100000 -s 10000 > window.100kb-10kb.chrZ-NC_053488.1.bed
echo -e 'NC_053453.1\t76187387' | bedtools makewindows -g - -w 50000 -s 5000 > window.50kb-5kb.chr1A-NC_053453.1.bed
echo -e 'NC_053450.1\t156035725' | bedtools makewindows -g - -w 50000 -s 5000 > window.50kb-5kb.chr2-NC_053450.1.bed
echo -e 'NC_053488.1\t90132487' | bedtools makewindows -g - -w 50000 -s 5000 > window.50kb-5kb.chrZ-NC_053488.1.bed
echo -e 'NC_053453.1\t76187387' | bedtools makewindows -g - -w 10000 -s 1000 > window.10kb-1kb.chr1A-NC_053453.1.bed
echo -e 'NC_053450.1\t156035725' | bedtools makewindows -g - -w 10000 -s 1000 > window.10kb-1kb.chr2-NC_053450.1.bed
echo -e 'NC_053488.1\t90132487' | bedtools makewindows -g - -w 10000 -s 1000 > window.10kb-1kb.chrZ-NC_053488.1.bed
```

2. Run sliding window analyses on each chromosome.
```
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr1A.vcf.gz --populations popmap.pixy --bed_file window.100kb-10kb.chr1A-NC_053453.1.bed --output_folder results_chrom-specific --output_prefix pixy_chr1A-NC_053453.1_100kb-10kb > pixy_chr1A-NC_053453.1_100kb-10kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr2.vcf.gz --populations popmap.pixy --bed_file window.100kb-10kb.chr2-NC_053450.1.bed --output_folder results_chrom-specific --output_prefix pixy_chr2-NC_053450.1_100kb-10kb > pixy_chr2-NC_053450.1_100kb-10kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --populations popmap.pixy --bed_file window.100kb-10kb.chrZ-NC_053488.1.bed --output_folder results_chrom-specific --output_prefix pixy_chrZ-NC_053488.1_100kb-10kb > pixy_chrZ-NC_053488.1_100kb-10kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr1A.vcf.gz --populations popmap.pixy --bed_file window.50kb-5kb.chr1A-NC_053453.1.bed --output_folder results_chrom-specific --output_prefix pixy_chr1A-NC_053453.1_50kb-5kb > pixy_chr1A-NC_053453.1_50kb-5kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr2.vcf.gz --populations popmap.pixy --bed_file window.50kb-5kb.chr2-NC_053450.1.bed --output_folder results_chrom-specific --output_prefix pixy_chr2-NC_053450.1_50kb-5kb > pixy_chr2-NC_053450.1_50kb-5kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --populations popmap.pixy --bed_file window.50kb-5kb.chrZ-NC_053488.1.bed --output_folder results_chrom-specific --output_prefix pixy_chrZ-NC_053488.1_50kb-5kb > pixy_chrZ-NC_053488.1_50kb-5kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr1A.vcf.gz --populations popmap.pixy --bed_file window.10kb-1kb.chr1A-NC_053453.1.bed --output_folder results_chrom-specific --output_prefix pixy_chr1A-NC_053453.1_10kb-1kb > pixy_chr1A-NC_053453.1_10kb-1kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr2.vcf.gz --populations popmap.pixy --bed_file window.10kb-1kb.chr2-NC_053450.1.bed --output_folder results_chrom-specific --output_prefix pixy_chr2-NC_053450.1_10kb-1kb > pixy_chr2-NC_053450.1_10kb-1kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --populations popmap.pixy --bed_file window.10kb-1kb.chrZ-NC_053488.1.bed --output_folder results_chrom-specific --output_prefix pixy_chrZ-NC_053488.1_10kb-1kb > pixy_chrZ-NC_053488.1_10kb-1kb.log
```

3. Run sliding window analyses on each chromosome between subspecies.
```
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr1A.vcf.gz --populations popmap.pixy.subspecies --bed_file window.100kb-10kb.chr1A-NC_053453.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chr1A-NC_053453.1_100kb-10kb > pixy_subspecies_chr1A-NC_053453.1_100kb-10kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr2.vcf.gz --populations popmap.pixy.subspecies --bed_file window.100kb-10kb.chr2-NC_053450.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chr2-NC_053450.1_100kb-10kb > pixy_subspecies_chr2-NC_053450.1_100kb-10kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --populations popmap.pixy.subspecies --bed_file window.100kb-10kb.chrZ-NC_053488.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chrZ-NC_053488.1_100kb-10kb > pixy_subspecies_chrZ-NC_053488.1_100kb-10kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr1A.vcf.gz --populations popmap.pixy.subspecies --bed_file window.50kb-5kb.chr1A-NC_053453.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chr1A-NC_053453.1_50kb-5kb > pixy_subspecies_chr1A-NC_053453.1_50kb-5kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr2.vcf.gz --populations popmap.pixy.subspecies --bed_file window.50kb-5kb.chr2-NC_053450.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chr2-NC_053450.1_50kb-5kb > pixy_subspecies_chr2-NC_053450.1_50kb-5kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --populations popmap.pixy.subspecies --bed_file window.50kb-5kb.chrZ-NC_053488.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chrZ-NC_053488.1_50kb-5kb > pixy_subspecies_chrZ-NC_053488.1_50kb-5kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr1A.vcf.gz --populations popmap.pixy.subspecies --bed_file window.10kb-1kb.chr1A-NC_053453.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chr1A-NC_053453.1_10kb-1kb > pixy_subspecies_chr1A-NC_053453.1_10kb-1kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chr2.vcf.gz --populations popmap.pixy.subspecies --bed_file window.10kb-1kb.chr2-NC_053450.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chr2-NC_053450.1_10kb-1kb > pixy_subspecies_chr2-NC_053450.1_10kb-1kb.log
pixy --n_cores 6 --stats pi dxy fst --vcf ~/hirundo_speciation_genomics/vcf/chrom-specific/hirundo_rustica+smithii.allsites.final.chrZ.vcf.gz --populations popmap.pixy.subspecies --bed_file window.10kb-1kb.chrZ-NC_053488.1.bed --output_folder results_chrom-specific --output_prefix pixy_subspecies_chrZ-NC_053488.1_10kb-1kb > pixy_subspecies_chrZ-NC_053488.1_10kb-1kb.log
```

### Format results

We'll concatenate windowed outputs for all of the chromosomes, which can be parsed further by population.

#### Format scaffold lists

1. Format complete scaffold list.
```
cat scaffold.auto+chrZ.list scaffold.chrW.list > scaffold.all.list
```
2. Format 'ordered' scaffold list (no unplaced scaffolds).
```
grep -v '-un' scaffold.all.list > scaffold.order.list
```
3. Run `concatenatePixy.sh` to combine results at different window resolutions.
```
sh concatenatePixy.sh
```
4. Run `concatenatePixyOrder.sh` to combine results for ordered chromosomes.
```
sh concatenatePixyOrder.sh
```
5. Run `concatenatePixyOrder_subspecies.sh`.
```
sh concatenatePixyOrder_subspecies.sh
```

### Fst values of significant GWA SNPS

To say whether regions of the genome strongly associated with traits have experienced divergent selection in the parental populations, we first need to explore the overlap between association peaks and Fst values by answering whether significant SNPs have higher Fst than non-significant SNPs and/or genome background Fst distributions.

#### Set up environment
```
cd ./analysis/gemma
mkdir fst
```

#### Convert Fst results to BED input format
```
awk 'BEGIN{OFS="\t"}{if ($1=="RU" && $2=="TY") print $3,$4-1,$5,$6}' ../pixy/pixy.all.order.fst.10kb.txt > ./fst/pixy.ruty.order.fst.10kb.bed
awk 'BEGIN{OFS="\t"}{if ($1=="RU" && $2=="GU") print $3,$4-1,$5,$6}' ../pixy/pixy.all.order.fst.10kb.txt > ./fst/pixy.rugu.order.fst.10kb.bed
awk 'BEGIN{OFS="\t"}{if ($1=="GU" && $2=="TY") print $3,$4-1,$5,$6}' ../pixy/pixy.all.order.fst.10kb.txt > ./fst/pixy.guty.order.fst.10kb.bed
```

#### Run bedtools intersect to output Fst values for significant GWA SNPs
```
for pop in ruty rugu guty; do echo -e "chrom\tsnp-start\tsnp-end\twald-p\twindow-start\twindow-end\tfst" > ./fst/gwas_full_lmm.full-breast-bright.assoc.sig.fst-${pop}.txt; bedtools intersect -a ./significant/gwas_full_lmm.full-breast-bright.assoc.sig.bed -b ./fst/pixy.$pop.order.fst.10kb.bed -wao | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7,$8}' | awk -F"\t" '!seen[$5, $6, $7]++' - >> ./fst/gwas_full_lmm.full-breast-bright.assoc.sig.fst-${pop}.txt; done
for pop in ruty rugu guty; do echo -e "chrom\tsnp-start\tsnp-end\twald-p\twindow-start\twindow-end\tfst" > ./fst/gwas_full_lmm.full-tail-streamer.assoc.sig.fst-${pop}.txt; bedtools intersect -a ./significant/gwas_full_lmm.full-tail-streamer.assoc.sig.bed -b ./fst/pixy.$pop.order.fst.10kb.bed -wao | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7,$8}' | awk -F"\t" '!seen[$5, $6, $7]++' - >> ./fst/gwas_full_lmm.full-tail-streamer.assoc.sig.fst-${pop}.txt; done
```

#### Run bedtools to output Fst values for non-significant GWA SNPs
```
for pop in ruty rugu guty; do echo -e "chrom\tsnp-start\tsnp-end\twald-p\twindow-start\twindow-end\tfst" > ./fst/gwas_full_lmm.full-breast-bright.assoc.nosig.fst-${pop}.txt; awk 'BEGIN{OFS="\t"}{if (-log($13)/log(10)<5) print $1,$3-1,$3,$13}' ./output/gwas_full_lmm.full-breast-bright.assoc.txt | grep -v 'NW_' | bedtools intersect -a - -b ./fst/pixy.$pop.order.fst.10kb.bed -wao | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7,$8}' | awk -F"\t" '!seen[$5, $6, $7]++' - >> ./fst/gwas_full_lmm.full-breast-bright.assoc.nosig.fst-${pop}.txt; done
for pop in ruty rugu guty; do echo -e "chrom\tsnp-start\tsnp-end\twald-p\twindow-start\twindow-end\tfst" > ./fst/gwas_full_lmm.full-tail-streamer.assoc.nosig.fst-${pop}.txt; awk 'BEGIN{OFS="\t"}{if (-log($13)/log(10)<5) print $1,$3-1,$3,$13}' ./output/gwas_full_lmm.full-tail-streamer.assoc.txt | grep -v 'NW_' | bedtools intersect -a - -b ./fst/pixy.$pop.order.fst.10kb.bed -wao | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7,$8}' | awk -F"\t" '!seen[$5, $6, $7]++' - >> ./fst/gwas_full_lmm.full-tail-streamer.assoc.nosig.fst-${pop}.txt; done
```

### Statistical analysis and comparison of Fst distributions

Run `./R/pixy.R` to summarize, plot, and statistically compare π, Fst, dxy, and recombination rate.

Run `./R/candidate_plotting_popgen_stats.R` to visualize population genetic summary statistics in candidate trait loci, along with additional statistics to detect signatures of selection (see below).

Run code blocks at the end of `./R/gemma_LMM.R` to summarize, plot, and statistically compare Fst distributions between significant and non-significant GWA SNPs.

[Back to top](#contents)


## Population branch statistics

It may be useful to further clarify which population(s) have been the targets of selection using population branch statistics to examine cases when one population has a much longer locus/region-specific branch length than the others, following [Yi et al. 2010](https://www.science.org/doi/full/10.1126/science.1190371): PBS = (T12 + T13 - T23)/2, where T = -log(1-Fst) for a given pair of populations.

We have Fst data required to calculate PBS in sliding windows from `pixy`. Run `./R/pbs.R` to perform calculations (also in `./R/candidate_plotting_popgen_stats.R`).

[Back to top](#contents)


## Tajima's D

We'll look for fluctuations in the allele frequency spectrum using Tajima's D, comparing Tajima and Watterson's θ. We'll estimate Tajima's D using `VCF-kit`, which has the ability to perform sliding window scans of the statistic. We do not want to limit analyses to SNPs with minor allele frequency cutoffs, since this will bias Watterson's estimator.

### Install VCF-kit

#### Set up environment
```
cd ~/hirundo_speciation_genomics/tmp/
mkdir vcf-kit-install
cd vcf-kit-install
```
#### Create and activate VCF-kit virtual environment
```
virtualenv -p /usr/bin/python3.8 vcf-kit-env
source ./vcf-kit-env/bin/activate
```
#### Install
```
pip install VCF-kit
```

### Set up input VCF data

#### Set up environment
```
cd ./analysis/
mkdir tajima
cd tajima
mkdir vcf
mkdir results
mkdir log
```

#### Format popmaps and chromosome list
```
popmap.rustica
popmap.tytleri
popmap.gutturalis
```

`chrom.list`

#### Extract VCFs for parental populations
```
bcftools view --threads 8 -S popmap.rustica -O z -o ./vcf/hirundo_rustica.parental.rustica.snps.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.vcf.gz
bcftools view --threads 8 -S popmap.tytleri -O z -o ./vcf/hirundo_rustica.parental.tytleri.snps.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.vcf.gz
bcftools view --threads 8 -S popmap.gutturalis -O z -o ./vcf/hirundo_rustica.parental.gutturalis.snps.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.vcf.gz
bcftools view --threads 8 -S popmap.rustica -O z -o ./vcf/hirundo_rustica.parental.rustica.snps.chrZ.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.vcf.gz
bcftools view --threads 8 -S popmap.tytleri -O z -o ./vcf/hirundo_rustica.parental.tytleri.snps.chrZ.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.vcf.gz
bcftools view --threads 8 -S popmap.gutturalis -O z -o ./vcf/hirundo_rustica.parental.gutturalis.snps.chrZ.vcf.gz ~/hirundo_speciation_genomics/vcf/hirundo_rustica+smithii.allsites.final.chrZ.snps.miss02.vcf.gz
```

### Calculate Tajima's D

#### Genome-wide analysis

Run `runTajima.sh` to calculate Tajima's D genome-wide.
```
source ~/hirundo_speciation_genomics/tmp/vcf-kit-install/vcf-kit-env/bin/activate
sh runTajima.sh
```

#### Focused analyses on specific chromosomes in sliding windows
```
bcftools view --threads 16 -r NC_053453.1 -O z ./vcf/hirundo_rustica.parental.rustica.snps.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.rustica.chr1A.50kb-5kb.txt
bcftools view --threads 16 -r NC_053453.1 -O z ./vcf/hirundo_rustica.parental.tytleri.snps.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.tytleri.chr1A.50kb-5kb.txt
bcftools view --threads 16 -r NC_053453.1 -O z ./vcf/hirundo_rustica.parental.gutturalis.snps.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.gutturalis.chr1A.50kb-5kb.txt
bcftools view --threads 16 -r NC_053488.1 -O z ./vcf/hirundo_rustica.parental.rustica.snps.chrZ.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.rustica.chrZ.50kb-5kb.txt
bcftools view --threads 16 -r NC_053488.1 -O z ./vcf/hirundo_rustica.parental.tytleri.snps.chrZ.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.tytleri.chrZ.50kb-5kb.txt
bcftools view --threads 16 -r NC_053488.1 -O z ./vcf/hirundo_rustica.parental.gutturalis.snps.chrZ.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.gutturalis.chrZ.50kb-5kb.txt
bcftools view --threads 16 -r NC_053450.1 -O z ./vcf/hirundo_rustica.parental.rustica.snps.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.rustica.chr2.50kb-5kb.txt
bcftools view --threads 16 -r NC_053450.1 -O z ./vcf/hirundo_rustica.parental.tytleri.snps.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.tytleri.chr2.50kb-5kb.txt
bcftools view --threads 16 -r NC_053450.1 -O z ./vcf/hirundo_rustica.parental.gutturalis.snps.vcf.gz | vk tajima 50000 5000 -  > ./results/tajima.gutturalis.chr2.50kb-5kb.txt
```

[Back to top](#contents)


## Haplotype statistics

We'll use haplotype diversity statistics to identify regions with patterns consistent with selection. First, we will phase variants for the parental populations. We'll then use the R package `rehh` to quantify haplotype statistics.

### Set up environment

```
cd ./analysis/
mkdir rehh
cd rehh
mkdir shapeit
cd shapeit
mkdir input
mkdir pirs
mkdir log
mkdir results
```

### Haplotype phasing

We'll perform statistical read-backed phasing using `SHAPEIT2`.

#### Retrieve SHAPEIT and extractPIRs executables

`wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz`
`tar -xf shapeit.v2.r904.glibcv2.17.linux.tar.gz`
`rm shapeit.v2.r904.glibcv2.17.linux.tar.gz`
`mv shapeit.v2.904.3.10.0-693.11.6.el7.x86_64 shapeit2`

`wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/files/extractPIRs.v1.r68.x86_64.tgz`
`tar -xf extractPIRs.v1.r68.x86_64.tgz`
`rm extractPIRs.v1.r68.x86_64.tgz`
`mv extractPIRs.v1.r68.x86_64 extractPIRs`

#### Prepare input data for phasing

1. Format population maps.
```
popmap.rustica
popmap.tytleri
popmap.gutturalis
```
2. Concatenate popmaps.
```
cat popmap.rustica popmap.tytleri popmap.gutturalis > popmap.all
```
3. Format `chromosome-scaffold.ordered.table.txt` and `chromosome-scaffold.table.txt`.
4. Run `parseParentalVCF.sh` to extract chromosome-specific VCFs for each parental population.
```
sh parseParentalVCF.sh chromosome-scaffold.ordered.table.txt
gunzip ./input/*.vcf.gz
```
5. Format bam list input files for `extractPIRs`.
```
sh makeBamlist.sh
rm bamlist.chr*-un*
```

#### Extract phase-informative reads

1. Format `chrom.list`.
2. Run `runExtractPIRs.sh` to extract phase informative reads.
```
sh runExtractPIRs.sh chrom.list
```

#### Assemble haplotypes using SHAPEIT2

We'll run `SHAPEIT2` using the settings:
* states = 1000
* burn = 200
* prune = 210
* main = 1000
* force

1. Run `runShapeitAssemble.sh` to assemble haplotypes.
```
sh runShapeitAssemble.sh chrom.list
```
2. Run `runShapeitConvert.sh` to convert haplotype output to VCF.
```
sh runShapeitConvert.sh chrom.list
```
3. Run `runCompressIndexVCFs.sh` to compress and index results.
```
sh runCompressIndexVCFs.sh chrom.list
```

#### Parse phased VCFs by population

Run `parsePhasedVCFs.sh` to extract parental population-specific VCFs.
```
sh parsePhasedVCFs.sh chrom.list
```

### Haplotype scans

We now have input data that we can analyze using `rehh`.
```
cd ./analysis/rehh/
```

Run `runScans.sh` to call the R script `rehhScans.R` on each chromosome, perform haplotype scans, and calculate iHS and xp-EHH.
```
sh runScans.sh
```

Concatenate results.
```
for pop in rustica tytleri gutturalis; do head -n1 ./results/iHS_chr1_${pop}.txt > ./results/iHS_all_${pop}.txt; for chrom in `cat chrom.list`; do tail -n+2 ./results/iHS_${chrom}_${pop}.txt >> ./results/iHS_all_${pop}.txt; done; done
for pop in rustica-tytleri rustica-gutturalis tytleri-gutturalis; do head -n1 ./results/XP-EHH_chr1_${pop}.txt > ./results/XP-EHH_all_${pop}.txt; for chrom in `cat chrom.list`; do tail -n+2 ./results/XP-EHH_${chrom}_${pop}.txt >> ./results/XP-EHH_all_${pop}.txt; done; done
```

[Back to top](#contents)


## Geographic cline analysis


## Genomic cline analysis


## Linkage disequilibrium: tests of genetic coupling


## Linkage disequilibrium: decay



















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





































