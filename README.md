# Sexual selection promotes reproductive isolation in barn swallows: Genomic analysis workflow

![Sexual selection promotes reproductive isolation in barn swallows](/cover-image.png "cover image")

This repository contains details on the data processing and analysis steps used to study the role of sexual selection in speciation through 1) the characterization of the genetic architecture of sexual signal traits in barn swallows (_Hirundo rustica_), 2) tests of divergent sexual selection in geographic isolation, 3) tests of selection against gene flow at trait loci in secondary contact, and 4) tests that selection has maintained associations between trait loci, known as genetic coupling, leading to reproductive isolation due to sexual selection. This workflow is a companion to the methods described in Schield et al. (_in review_).

Lists and reference files can be found in various analysis directories, as can Shell and Python scripts used in elements of the analysis workflow. R scripts are in the `R` directory. Note that you may need to adjust the organization of path/file locations to suit your environment. This workflow assumes that you return to the main working directory `~/hirundo_speciation_genomics` after completing major analysis sections. Adjust the path before the main working directory as necessary. Walkthrough of software installations will take place in `~/hirundo_speciation_genomics/tmp/`.

Feel free to contact me at drew.schield[at]colorado.edu with any questions.

## Contents

* [Software and dependencies](#software-and-dependencies)
* [Genome data processing and variant calling](#genome-data-processing-and-variant-calling)
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
* [bcftools](http://www.htslib.org/)
* [Plink](https://www.cog-genomics.org/plink/)
* [vcftools](https://vcftools.github.io/index.html)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [ADMIXTURE](https://dalexander.github.io/admixture/publications.html)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [rehh](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html)
* [MashMap](https://github.com/marbl/MashMap)
* [Repeatmasker](https://www.repeatmasker.org/)
* [R](https://cran.r-project.org/)

Note, I installed a number of these programs to my [conda](https://docs.conda.io/en/latest/) environment, or installed via a virtual environment (details below).

[Back to top](#contents)

## Genome data processing and variant calling

In progress...

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










































