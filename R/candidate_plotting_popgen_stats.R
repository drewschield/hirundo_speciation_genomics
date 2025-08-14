############################################################################
# Plotting outlier regions associated with phenotypes - selection statistics
############################################################################

# We identified genotype X phenotype associations using GEMMA and estimated
# π, dxy, and Fst within and between populations and their respective
# hybrid zones using Pixy.

# Here, we will focus on plotting results for genomic regions housing
# genes intersecting or very near (i.e., +/- 20kb) strongly signficant SNPs
# from analysis in GEMMA. A number of genes associated with breast
# brightness, for example, are known candidate genes involved in
# melanogenesis and pigmentation.These scans will include pi, dxy, Fst, PBS,
# Tajima's D, and iHS/xpEHH.

# This script contains commands for parsing and plotting data from 
# population genetic summary statistics in these focal genomic regions.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/candidate_plotting/')

library(data.table)
library(tidyverse)
library(scales)
library(ggplot2)

### General reference-------------------------------------------------------

# Most strong association peaks with phenotypes are concentrated on
# chromosomes 1A, 2, and the Z chromosome. These are the chromosome-
# scaffold associations:

# Chromosome  Scaffold
# Chr. 1A     NC_053453.1
# Chr. 2      NC_053450.1
# Z chr.      NC_053488.1

### General reference-------------------------------------------------------

# Strong association peaks with phenotypes are concentrated on
# chromosomes 1A, 2, the Z chromosome. These are the chromosome-
# scaffold associations:

# Chromosome  Scaffold
# Chr. 1A     NC_053453.1
# Chr. 2      NC_053450.1
# Z chr.      NC_053488.1

### Read in candidate gene data----------------------------------------------

# Read data
cand.vc <- read.table('./coordinates_candidate_gene_breast-brightness.txt',header=T)
cand.ts <- read.table('./coordinates_candidate_gene_tail-streamer.txt',header=T)

# Parse strand
cand.vc.for <- cand.vc[which(cand.vc$strand=='+'),]
cand.vc.rev <- cand.vc[which(cand.vc$strand=='-'),]

cand.ts.for <- cand.ts[which(cand.ts$strand=='+'),]
cand.ts.rev <- cand.ts[which(cand.ts$strand=='-'),]

### Read in general gene/exon data-=------------------------------------------

# Read data
genes <- read.table('./GCF_015227805.1_bHirRus1.pri.v2_genomic.gene.sort.txt',header=T)
exons <- read.table('./GCF_015227805.1_bHirRus1.pri.v2_genomic.exon.sort.txt',header=T)

# Parse strand
genes.for <- genes[which(genes$strand=='+'),]
genes.rev <- genes[which(genes$strand=='-'),]

exons.for <- exons[which(exons$strand=='+'),]
exons.rev <- exons[which(exons$strand=='-'),]

### Read in population genetic statistics (smoothed sliding windows)----------

# Fst
## Read data
fst.100kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_100kb-10kb_fst.txt',header=T)
fst.100kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_100kb-10kb_fst.txt',header=T)
fst.100kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_100kb-10kb_fst.txt',header=T)
fst.50kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_50kb-5kb_fst.txt',header=T)
fst.50kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_50kb-5kb_fst.txt',header=T)
fst.50kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_50kb-5kb_fst.txt',header=T)
fst.10kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_10kb-1kb_fst.txt',header=T)
fst.10kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_10kb-1kb_fst.txt',header=T)
fst.10kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_10kb-1kb_fst.txt',header=T)
## Parse populations
fst.100kb.1a.ruty <- fst.100kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.100kb.1a.rugu <- fst.100kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.100kb.1a.guty <- fst.100kb.1a %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.100kb.2.ruty <- fst.100kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.100kb.2.rugu <- fst.100kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.100kb.2.guty <- fst.100kb.2 %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.100kb.z.ruty <- fst.100kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.100kb.z.rugu <- fst.100kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.100kb.z.guty <- fst.100kb.z %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.50kb.1a.ruty <- fst.50kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.50kb.1a.rugu <- fst.50kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.50kb.1a.guty <- fst.50kb.1a %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.50kb.2.ruty <- fst.50kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.50kb.2.rugu <- fst.50kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.50kb.2.guty <- fst.50kb.2 %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.50kb.z.ruty <- fst.50kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.50kb.z.rugu <- fst.50kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.50kb.z.guty <- fst.50kb.z %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.10kb.1a.ruty <- fst.10kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.10kb.1a.rugu <- fst.10kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.10kb.1a.guty <- fst.10kb.1a %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.10kb.2.ruty <- fst.10kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.10kb.2.rugu <- fst.10kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.10kb.2.guty <- fst.10kb.2 %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.10kb.z.ruty <- fst.10kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.10kb.z.rugu <- fst.10kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
fst.10kb.z.guty <- fst.10kb.z %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))

# dxy
## read data
dxy.100kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_100kb-10kb_dxy.txt',header=T)
dxy.100kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_100kb-10kb_dxy.txt',header=T)
dxy.100kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_100kb-10kb_dxy.txt',header=T)
dxy.50kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_50kb-5kb_dxy.txt',header=T)
dxy.50kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_50kb-5kb_dxy.txt',header=T)
dxy.50kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_50kb-5kb_dxy.txt',header=T)
dxy.10kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_10kb-1kb_dxy.txt',header=T)
dxy.10kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_10kb-1kb_dxy.txt',header=T)
dxy.10kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_10kb-1kb_dxy.txt',header=T)
## parse populations
dxy.100kb.1a.ruty <- dxy.100kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.100kb.1a.rugu <- dxy.100kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.100kb.1a.guty <- dxy.100kb.1a %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.100kb.2.ruty <- dxy.100kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.100kb.2.rugu <- dxy.100kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.100kb.2.guty <- dxy.100kb.2 %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.100kb.z.ruty <- dxy.100kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.100kb.z.rugu <- dxy.100kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.100kb.z.guty <- dxy.100kb.z %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.50kb.1a.ruty <- dxy.50kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.50kb.1a.rugu <- dxy.50kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.50kb.1a.guty <- dxy.50kb.1a %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.50kb.2.ruty <- dxy.50kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.50kb.2.rugu <- dxy.50kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.50kb.2.guty <- dxy.50kb.2 %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.50kb.z.ruty <- dxy.50kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.50kb.z.rugu <- dxy.50kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.50kb.z.guty <- dxy.50kb.z %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.10kb.1a.ruty <- dxy.10kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.10kb.1a.rugu <- dxy.10kb.1a %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.10kb.1a.guty <- dxy.10kb.1a %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.10kb.2.ruty <- dxy.10kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.10kb.2.rugu <- dxy.10kb.2 %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.10kb.2.guty <- dxy.10kb.2 %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
dxy.10kb.z.ruty <- dxy.10kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.10kb.z.rugu <- dxy.10kb.z %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'GU'))
dxy.10kb.z.guty <- dxy.10kb.z %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))

# pi
## read data
pi.100kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_100kb-10kb_pi.txt',header=T)
pi.100kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_100kb-10kb_pi.txt',header=T)
pi.100kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_100kb-10kb_pi.txt',header=T)
pi.50kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_50kb-5kb_pi.txt',header=T)
pi.50kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_50kb-5kb_pi.txt',header=T)
pi.50kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_50kb-5kb_pi.txt',header=T)
pi.10kb.1a <- read.table('../pixy/results_chrom-specific/pixy_chr1A-NC_053453.1_10kb-1kb_pi.txt',header=T)
pi.10kb.2 <- read.table('../pixy/results_chrom-specific/pixy_chr2-NC_053450.1_10kb-1kb_pi.txt',header=T)
pi.10kb.z <- read.table('../pixy/results_chrom-specific/pixy_chrZ-NC_053488.1_10kb-1kb_pi.txt',header=T)
## parse populations
pi.100kb.1a.ru <- pi.100kb.1a %>% filter(str_detect(pop, 'RU'))
pi.100kb.1a.ty <- pi.100kb.1a %>% filter(str_detect(pop, 'TY'))
pi.100kb.1a.gu <- pi.100kb.1a %>% filter(str_detect(pop, 'GU'))
pi.100kb.2.ru <- pi.100kb.2 %>% filter(str_detect(pop, 'RU'))
pi.100kb.2.ty <- pi.100kb.2 %>% filter(str_detect(pop, 'TY'))
pi.100kb.2.gu <- pi.100kb.2 %>% filter(str_detect(pop, 'GU'))
pi.100kb.z.ru <- pi.100kb.z %>% filter(str_detect(pop, 'RU'))
pi.100kb.z.ty <- pi.100kb.z %>% filter(str_detect(pop, 'TY'))
pi.100kb.z.gu <- pi.100kb.z %>% filter(str_detect(pop, 'GU'))
pi.50kb.1a.ru <- pi.50kb.1a %>% filter(str_detect(pop, 'RU'))
pi.50kb.1a.ty <- pi.50kb.1a %>% filter(str_detect(pop, 'TY'))
pi.50kb.1a.gu <- pi.50kb.1a %>% filter(str_detect(pop, 'GU'))
pi.50kb.2.ru <- pi.50kb.2 %>% filter(str_detect(pop, 'RU'))
pi.50kb.2.ty <- pi.50kb.2 %>% filter(str_detect(pop, 'TY'))
pi.50kb.2.gu <- pi.50kb.2 %>% filter(str_detect(pop, 'GU'))
pi.50kb.z.ru <- pi.50kb.z %>% filter(str_detect(pop, 'RU'))
pi.50kb.z.ty <- pi.50kb.z %>% filter(str_detect(pop, 'TY'))
pi.50kb.z.gu <- pi.50kb.z %>% filter(str_detect(pop, 'GU'))
pi.10kb.1a.ru <- pi.10kb.1a %>% filter(str_detect(pop, 'RU'))
pi.10kb.1a.ty <- pi.10kb.1a %>% filter(str_detect(pop, 'TY'))
pi.10kb.1a.gu <- pi.10kb.1a %>% filter(str_detect(pop, 'GU'))
pi.10kb.2.ru <- pi.10kb.2 %>% filter(str_detect(pop, 'RU'))
pi.10kb.2.ty <- pi.10kb.2 %>% filter(str_detect(pop, 'TY'))
pi.10kb.2.gu <- pi.10kb.2 %>% filter(str_detect(pop, 'GU'))
pi.10kb.z.ru <- pi.10kb.z %>% filter(str_detect(pop, 'RU'))
pi.10kb.z.ty <- pi.10kb.z %>% filter(str_detect(pop, 'TY'))
pi.10kb.z.gu <- pi.10kb.z %>% filter(str_detect(pop, 'GU'))

# PBS
# Note - full details of calculations are in 'PBS.R'.
pbs.100kb <- read.table('../pbs/pbs.100kb.txt',header=T)
pbs.10kb <- read.table('../pbs/pbs.10kb.txt',header=T)

pbs.100kb.1a <- read.table('../pbs/pbs.100kb.chr1A.txt',header=T)
pbs.50kb.1a <- read.table('../pbs/pbs.50kb.chr1A.txt',header=T)
pbs.10kb.1a <- read.table('../pbs/pbs.10kb.chr1A.txt',header=T)

pbs.100kb.2 <- read.table('../pbs/pbs.100kb.chr2.txt',header=T)
pbs.50kb.2 <- read.table('../pbs/pbs.50kb.chr2.txt',header=T)
pbs.10kb.2 <- read.table('../pbs/pbs.10kb.chr2.txt',header=T)

pbs.100kb.z <- read.table('../pbs/pbs.100kb.chrZ.txt',header=T)
pbs.50kb.z <- read.table('../pbs/pbs.50kb.chrZ.txt',header=T)
pbs.10kb.z <- read.table('../pbs/pbs.10kb.chrZ.txt',header=T)

# Tajima's D
taj.10kb.ru <- read.table('../tajima/tajima.rustica.10kb.txt',header=T)
taj.10kb.ty <- read.table('../tajima/tajima.tytleri.10kb.txt',header=T)
taj.10kb.gu <- read.table('../tajima/tajima.gutturalis.10kb.txt',header=T)

taj.50kb.1a.ru <- read.table('../tajima/tajima.rustica.chr1A.50kb-5kb.txt',header=T)
taj.50kb.1a.ty <- read.table('../tajima/tajima.tytleri.chr1A.50kb-5kb.txt',header=T)
taj.50kb.1a.gu <- read.table('../tajima/tajima.gutturalis.chr1A.50kb-5kb.txt',header=T)
taj.50kb.2.ru <- read.table('../tajima/tajima.rustica.chr2.50kb-5kb.txt',header=T)
taj.50kb.2.ty <- read.table('../tajima/tajima.tytleri.chr2.50kb-5kb.txt',header=T)
taj.50kb.2.gu <- read.table('../tajima/tajima.gutturalis.chr2.50kb-5kb.txt',header=T)
taj.50kb.z.ru <- read.table('../tajima/tajima.rustica.chrZ.50kb-5kb.txt',header=T)
taj.50kb.z.ty <- read.table('../tajima/tajima.tytleri.chrZ.50kb-5kb.txt',header=T)
taj.50kb.z.gu <- read.table('../tajima/tajima.gutturalis.chrZ.50kb-5kb.txt',header=T)

# iHS
ihs.all.ru <- read.table('../rehh/results/iHS_all_rustica.txt',header=T)
ihs.all.ty <- read.table('../rehh/results/iHS_all_tytleri.txt',header=T)
ihs.all.gu <- read.table('../rehh/results/iHS_all_gutturalis.txt',header=T)

## Get genome-wide quantiles
quantile(ihs.all.ru$IHS,c(0.01,0.99),na.rm=T)
## -2.833746, 2.211135
quantile(ihs.all.ty$IHS,c(0.01,0.99),na.rm=T)
## -2.864729, 2.169999
quantile(ihs.all.gu$IHS,c(0.01,0.99),na.rm=T)
## -2.737471, 2.211431

ihs.1a.ru <- read.table('../rehh/results/iHS_chr1A_rustica.txt',header=T)
ihs.1a.ty <- read.table('../rehh/results/iHS_chr1A_tytleri.txt',header=T)
ihs.1a.gu <- read.table('../rehh/results/iHS_chr1A_gutturalis.txt',header=T)
ihs.2.ru <- read.table('../rehh/results/iHS_chr2_rustica.txt',header=T)
ihs.2.ty <- read.table('../rehh/results/iHS_chr2_tytleri.txt',header=T)
ihs.2.gu <- read.table('../rehh/results/iHS_chr2_gutturalis.txt',header=T)
ihs.z.ru <- read.table('../rehh/results/iHS_chrZ_rustica.txt',header=T)
ihs.z.ty <- read.table('../rehh/results/iHS_chrZ_tytleri.txt',header=T)
ihs.z.gu <- read.table('../rehh/results/iHS_chrZ_gutturalis.txt',header=T)

# XP-EHH
xpehh.all.ruty <- read.table('../rehh/results/XP-EHH_all_rustica-tytleri.txt',header=T)
xpehh.all.rugu <- read.table('../rehh/results/XP-EHH_all_rustica-gutturalis.txt',header=T)
xpehh.all.guty <- read.table('../rehh/results/XP-EHH_all_tytleri-gutturalis.txt',header=T)
## Get genome-wide quantiles
quantile(xpehh.all.ruty$XPEHH_rustica_tytleri,c(0.01,0.99),na.rm=T)
## -2.497729, 2.477378
quantile(xpehh.all.rugu$XPEHH_rustica_gutturalis,c(0.01,0.99),na.rm=T)
## -2.570383, 2.417505
quantile(xpehh.all.guty$XPEHH_tytleri_gutturalis,c(0.01,0.99),na.rm=T)
## -2.542992, 2.447439

xpehh.1a.ruty <- read.table('../rehh/results/XP-EHH_chr1A_rustica-tytleri.txt',header=T)
xpehh.1a.rugu <- read.table('../rehh/results/XP-EHH_chr1A_rustica-gutturalis.txt',header=T)
xpehh.1a.guty <- read.table('../rehh/results/XP-EHH_chr1A_tytleri-gutturalis.txt',header=T)
xpehh.2.ruty <- read.table('../rehh/results/XP-EHH_chr2_rustica-tytleri.txt',header=T)
xpehh.2.rugu <- read.table('../rehh/results/XP-EHH_chr2_rustica-gutturalis.txt',header=T)
xpehh.2.guty <- read.table('../rehh/results/XP-EHH_chr2_tytleri-gutturalis.txt',header=T)
xpehh.z.ruty <- read.table('../rehh/results/XP-EHH_chrZ_rustica-tytleri.txt',header=T)
xpehh.z.rugu <- read.table('../rehh/results/XP-EHH_chrZ_rustica-gutturalis.txt',header=T)
xpehh.z.guty <- read.table('../rehh/results/XP-EHH_chrZ_tytleri-gutturalis.txt',header=T)

xpehh.100kb.1a.ruty <- read.table('../rehh/results_windowed/XP-EHH_chr1A_rustica-tytleri.window.100kb-10kb.txt',header=T)
xpehh.100kb.1a.rugu <- read.table('../rehh/results_windowed/XP-EHH_chr1A_rustica-gutturalis.window.100kb-10kb.txt',header=T)
xpehh.100kb.1a.guty <- read.table('../rehh/results_windowed/XP-EHH_chr1A_tytleri-gutturalis.window.100kb-10kb.txt',header=T)
xpehh.100kb.2.ruty <- read.table('../rehh/results_windowed/XP-EHH_chr2_rustica-tytleri.window.100kb-10kb.txt',header=T)
xpehh.100kb.2.rugu <- read.table('../rehh/results_windowed/XP-EHH_chr2_rustica-gutturalis.window.100kb-10kb.txt',header=T)
xpehh.100kb.2.guty <- read.table('../rehh/results_windowed/XP-EHH_chr2_tytleri-gutturalis.window.100kb-10kb.txt',header=T)
xpehh.100kb.z.ruty <- read.table('../rehh/results_windowed/XP-EHH_chrZ_rustica-tytleri.window.100kb-10kb.txt',header=T)
xpehh.100kb.z.rugu <- read.table('../rehh/results_windowed/XP-EHH_chrZ_rustica-gutturalis.window.100kb-10kb.txt',header=T)
xpehh.100kb.z.guty <- read.table('../rehh/results_windowed/XP-EHH_chrZ_tytleri-gutturalis.window.100kb-10kb.txt',header=T)

xpehh.50kb.1a.ruty <- read.table('../rehh/results_windowed/XP-EHH_chr1A_rustica-tytleri.window.50kb-5kb.txt',header=T)
xpehh.50kb.1a.rugu <- read.table('../rehh/results_windowed/XP-EHH_chr1A_rustica-gutturalis.window.50kb-5kb.txt',header=T)
xpehh.50kb.1a.guty <- read.table('../rehh/results_windowed/XP-EHH_chr1A_tytleri-gutturalis.window.50kb-5kb.txt',header=T)
xpehh.50kb.2.ruty <- read.table('../rehh/results_windowed/XP-EHH_chr2_rustica-tytleri.window.50kb-5kb.txt',header=T)
xpehh.50kb.2.rugu <- read.table('../rehh/results_windowed/XP-EHH_chr2_rustica-gutturalis.window.50kb-5kb.txt',header=T)
xpehh.50kb.2.guty <- read.table('../rehh/results_windowed/XP-EHH_chr2_tytleri-gutturalis.window.50kb-5kb.txt',header=T)
xpehh.50kb.z.ruty <- read.table('../rehh/results_windowed/XP-EHH_chrZ_rustica-tytleri.window.50kb-5kb.txt',header=T)
xpehh.50kb.z.rugu <- read.table('../rehh/results_windowed/XP-EHH_chrZ_rustica-gutturalis.window.50kb-5kb.txt',header=T)
xpehh.50kb.z.guty <- read.table('../rehh/results_windowed/XP-EHH_chrZ_tytleri-gutturalis.window.50kb-5kb.txt',header=T)

# Recombination rate
rho.100kb.1a <- read.table('../pyrho/rustica.rmap.chr1A-NC_053453.1.100kb-10kb.txt',header=T)
rho.100kb.2 <- read.table('../pyrho/rustica.rmap.chr2-NC_053450.1.100kb-10kb.txt',header=T)
rho.100kb.z <- read.table('../pyrho/rustica.rmap.chrZ-NC_053488.1.100kb-10kb.txt',header=T)

rho.all <- read.table('../pyrho/rustica.rmap.100kb.txt',header=T)
rho.auto <- rho.all %>% filter(!str_detect(chrom,'NC_053488.1'))

## Do quick comparison of recombination rate in trait loci compared to autosome / Z chromosome backgrounds
rho.1a.foc <- rho.100kb.1a %>% filter(start >= 4.375e+07 & end <= 4.5e+07 | start >= 4.62e+07 & end <= 4.69e+07)
rho.2.foc <- rho.100kb.2 %>% filter(start >= 9.78e+07 & end <= 9.855e+07 | start >= 1.e+08 & end <= 1.015e+08)
rho.z.foc <- rho.100kb.z %>% filter(start >= 1.85e+07 & end <= 1.953e+07 | 
                                                                 start >= 3.275e+07 & end <= 3.375e+07 | 
                                                                 start >= 3.5e+07 & end <= 3.69e+07 | 
                                                                 start >= 3.5e+07 & end <= 3.69e+07 | 
                                                                 start >= 4.36e+07 & end <= 4.45e+07 | 
                                                                 start >= 4.46e+07 & end <= 4.83e+07)
rho.auto.foc <- rbind(rho.1a.foc,rho.2.foc)

mean(rho.auto.foc$rate)
mean(as.numeric(rho.auto$rate),na.rm=T)
wilcox.test(rho.auto.foc$rate,as.numeric(rho.auto$rate))

mean(rho.z.foc$rate)
mean(as.numeric(rho.100kb.z$rate),na.rm=T)
wilcox.test(rho.z.foc$rate,as.numeric(rho.100kb.z$rate))

### Plot chromosome-wide results----------------------------------------------

n <- 10
fst.1mb.1a.ruty <- aggregate(fst.100kb.1a.ruty, list(rep(1:(nrow(fst.100kb.1a.ruty) %/% n + 1), each = n, len = nrow(fst.100kb.1a.ruty))), mean)[-1];
dxy.1mb.1a.ruty <- aggregate(dxy.100kb.1a.ruty, list(rep(1:(nrow(dxy.100kb.1a.ruty) %/% n + 1), each = n, len = nrow(dxy.100kb.1a.ruty))), mean)[-1];
pi.1mb.1a.ru <- aggregate(pi.100kb.1a.ru, list(rep(1:(nrow(pi.100kb.1a.ru) %/% n + 1), each = n, len = nrow(pi.100kb.1a.ru))), mean)[-1];
rho.1mb.1a <- aggregate(rho.100kb.1a, list(rep(1:(nrow(rho.100kb.1a) %/% n + 1), each = n, len = nrow(rho.100kb.1a))), mean)[-1];

fst.1mb.z.ruty <- aggregate(fst.100kb.z.ruty, list(rep(1:(nrow(fst.100kb.z.ruty) %/% n + 1), each = n, len = nrow(fst.100kb.z.ruty))), mean)[-1];
dxy.1mb.z.ruty <- aggregate(dxy.100kb.z.ruty, list(rep(1:(nrow(dxy.100kb.z.ruty) %/% n + 1), each = n, len = nrow(dxy.100kb.z.ruty))), mean)[-1];
pi.1mb.z.ru <- aggregate(pi.100kb.z.ru, list(rep(1:(nrow(pi.100kb.z.ru) %/% n + 1), each = n, len = nrow(pi.100kb.z.ru))), mean)[-1];
rho.1mb.z <- aggregate(rho.100kb.z, list(rep(1:(nrow(rho.100kb.z) %/% n + 1), each = n, len = nrow(rho.100kb.z))), mean)[-1];

# Export at 10.5 x 7
par(mfrow=c(4,1))
plot(as.integer(row.names(fst.1mb.1a.ruty))*10,fst.1mb.1a.ruty$avg_wc_fst,type='l',lwd=2,ylim=c(0,1),xlim=c(0,156035725/10000),axes=T,ylab='Fst',xlab="")
for (gene in cand.br.1a){abline(v=(cand.br.1a$start.gene+(cand.br.1a$end.gene-cand.br.1a$start.gene)/2)/10000,lty=2)}
plot(as.integer(row.names(dxy.1mb.1a.ruty))*10,dxy.1mb.1a.ruty$avg_dxy,type='l',lwd=2,ylim=c(0,0.025),xlim=c(0,156035725/10000),axes=T,ylab='dxy',xlab="")
for (gene in cand.br.1a){abline(v=(cand.br.1a$start.gene+(cand.br.1a$end.gene-cand.br.1a$start.gene)/2)/10000,lty=2)}
plot(as.integer(row.names(pi.1mb.1a.ru))*10,pi.1mb.1a.ru$avg_pi,type='l',lwd=2,ylim=c(0,0.025),xlim=c(0,156035725/10000),axes=T,ylab='π',xlab="")
for (gene in cand.br.1a){abline(v=(cand.br.1a$start.gene+(cand.br.1a$end.gene-cand.br.1a$start.gene)/2)/10000,lty=2)}
plot(as.integer(row.names(rho.1mb.1a))*10,rho.1mb.1a$rate,type='l',lwd=2,ylim=c(0,2.5e-06),xlim=c(0,156035725/10000),axes=T,ylab='Recombination Rate',xlab="")
for (gene in cand.br.1a){abline(v=(cand.br.1a$start.gene+(cand.br.1a$end.gene-cand.br.1a$start.gene)/2)/10000,lty=2)}

par(mfrow=c(4,1))
plot(as.integer(row.names(fst.1mb.z.ruty))*10,fst.1mb.z.ruty$avg_wc_fst,type='l',lwd=2,ylim=c(0,1),xlim=c(0,156035725/10000),axes=T,ylab='Fst',xlab="")
for (gene in cand.br.z){abline(v=(cand.br.z$start.gene+(cand.br.z$end.gene-cand.br.z$start.gene)/2)/10000,lty=2)}
plot(as.integer(row.names(dxy.1mb.z.ruty))*10,dxy.1mb.z.ruty$avg_dxy,type='l',lwd=2,ylim=c(0,0.014),xlim=c(0,156035725/10000),axes=T,ylab='dxy',xlab="")
for (gene in cand.br.z){abline(v=(cand.br.z$start.gene+(cand.br.z$end.gene-cand.br.z$start.gene)/2)/10000,lty=2)}
plot(as.integer(row.names(pi.1mb.z.ru))*10,pi.1mb.z.ru$avg_pi,type='l',lwd=2,ylim=c(0,0.014),xlim=c(0,156035725/10000),axes=T,ylab='π',xlab="")
for (gene in cand.br.z){abline(v=(cand.br.z$start.gene+(cand.br.z$end.gene-cand.br.z$start.gene)/2)/10000,lty=2)}
plot(as.integer(row.names(rho.1mb.z))*10,rho.1mb.z$rate,type='l',lwd=2,ylim=c(0,10e-08),xlim=c(0,156035725/10000),axes=T,ylab='Recombination Rate',xlab="")
for (gene in cand.br.z){abline(v=(cand.br.z$start.gene+(cand.br.z$end.gene-cand.br.z$start.gene)/2)/10000,lty=2)}

### Plot Chromosome 1A demo (for main figure)---------------------------------

# Output at 10 x 4.5

par(mfrow=c(4,1))
plot(as.integer(row.names(fst.1mb.1a.ruty))*10,fst.1mb.1a.ruty$avg_wc_fst,type='l',lwd=2,ylim=c(0,0.4),axes=T,ylab='Fst',xlab="",col='#E28026')
plot(as.integer(row.names(dxy.1mb.1a.ruty))*10,dxy.1mb.1a.ruty$avg_dxy,type='l',lwd=2,ylim=c(0,0.012),axes=T,ylab='dxy',xlab="",col='#E28026')
plot(as.integer(row.names(pi.1mb.1a.ru))*10,pi.1mb.1a.ru$avg_pi,type='l',lwd=2,ylim=c(0,0.012),axes=T,ylab='π',xlab="",col='#B03160')
plot(as.integer(row.names(rho.1mb.1a))*10,rho.1mb.1a$rate,type='l',lwd=2,ylim=c(0,2.25e-06),axes=T,ylab='Recombination Rate',xlab="Chromosome Position (Mb)",col='#B03160')

### Plot KITLG region demo (for main figure)----------------------------------

# Plot PBS, π, Tajima's D, and xp-EHH (rustica vs tytleri) scans

## Set region & objects
x <- c(4.375e+07,4.5e+07)
y <- 0.5
y2 <- y-0.1
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053453.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053453.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.1a.ruty.foc <- xpehh.1a.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)

par(mfrow=c(4,1))
## PBS scan
y <- 0.55
y2 <- y-0.1
plot(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.ru,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.6))
lines(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
y <- 0.0027
y2 <- 0.0022
plot(pi.50kb.1a.ru$window_pos_1+(pi.50kb.1a.ru$window_pos_2 - pi.50kb.1a.ru$window_pos_1)/2,pi.50kb.1a.ru$avg_pi,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.0025))
lines(pi.50kb.1a.ty$window_pos_1+(pi.50kb.1a.ty$window_pos_2 - pi.50kb.1a.ty$window_pos_1)/2,pi.50kb.1a.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.1a.gu$window_pos_1+(pi.50kb.1a.gu$window_pos_2 - pi.50kb.1a.gu$window_pos_1)/2,pi.50kb.1a.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
y <- 4
y2 <- 3.5
plot(taj.50kb.1a.ru$BIN_START+(taj.50kb.1a.ru$BIN_END - taj.50kb.1a.ru$BIN_START)/2,taj.50kb.1a.ru$TajimaD,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.5,1))
lines(taj.50kb.1a.ty$BIN_START+(taj.50kb.1a.ty$BIN_END - taj.50kb.1a.ty$BIN_START)/2,taj.50kb.1a.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.1a.gu$BIN_START+(taj.50kb.1a.gu$BIN_END - taj.50kb.1a.gu$BIN_START)/2,taj.50kb.1a.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col='grey')
## xp-EHH scan
y <- 11
y2 <- 9.5
plot(xpehh.1a.ruty.foc$POSITION,xpehh.1a.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col='#B03160',xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,lwd=1,col='grey')


# Plot boxplots to match scans for KITLG

## Set random background points and focal data
pbs.all.ru.rand <- pbs.10kb[sample(nrow(pbs.10kb), 500), ]
pbs.1a.ru.back <- pbs.10kb %>% filter(str_detect(fst.10kb.ruty.chromosome, 'NC_053453.1'))
pbs.1a.ru.rand <- pbs.1a.ru.back[sample(nrow(pbs.1a.ru.back), 500),]
pbs.1a.ru.foc <- pbs.10kb %>% filter(str_detect(fst.10kb.ruty.chromosome, 'NC_053453.1') & (fst.10kb.ruty.window_pos_1 >= x[1] & fst.10kb.ruty.window_pos_2 <= x[2]))

pi.all.ru <- read.table('../pixy/pixy.all.order.pi.10kb.txt',header=T)
pi.all.ru <- pi.all.ru %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop,'RU'))
pi.all.ru.rand <- pi.all.ru[sample(nrow(pi.all.ru), 500),]
pi.1a.ru.back <- pi.all.ru %>% filter(str_detect(chromosome,'NC_053453.1'))
pi.1a.ru.rand <- pi.1a.ru.back[sample(nrow(pi.1a.ru.back), 500),]
pi.1a.ru.foc <- pi.all.ru %>% filter(str_detect(chromosome,'NC_053453.1') & (window_pos_1 >= x[1] & window_pos_2 <= x[2]))

taj.all.ru.rand <- taj.10kb.ru[sample(nrow(taj.10kb.ru),500),]
taj.1a.ru.back <- taj.10kb.ru %>% filter(str_detect(CHROM,'NC_053453.1'))
taj.1a.ru.rand <- taj.1a.ru.back[sample(nrow(taj.1a.ru.back),500),]
taj.1a.ru.foc <- taj.10kb.ru %>% filter(str_detect(CHROM,'NC_053453.1') & (BIN_START >= x[1] & BIN_END <= x[2]))

ggplot() +
  geom_boxplot(data=pbs.1a.ru.foc, aes(x = 1, y = pbs.10kb.ru),outlier.shape = NA) +
  geom_boxplot(data=pbs.1a.ru.back, aes(x = 2, y = pbs.10kb.ru),outlier.shape = NA) +
  geom_boxplot(data=pbs.10kb, aes(x = 3, y = pbs.10kb.ru),outlier.shape = NA) +
  geom_jitter(data=pbs.1a.ru.foc, aes(x = 1, y = pbs.10kb.ru),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=pbs.1a.ru.rand, aes(x = 2, y = pbs.10kb.ru),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=pbs.all.ru.rand, aes(x = 3, y = pbs.10kb.ru),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  ylim(0,0.6) +
  theme_classic()

ggplot() +
  geom_boxplot(data=pi.1a.ru.foc, aes(x = 1, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.1a.ru.back, aes(x = 2, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.all.ru, aes(x = 3, y = avg_pi),outlier.shape = NA) +
  geom_jitter(data=pi.1a.ru.foc, aes(x = 1, y = avg_pi),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=pi.1a.ru.rand, aes(x = 2, y = avg_pi),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=pi.all.ru.rand, aes(x = 3, y = avg_pi),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  ylim(0,0.01) +
  theme_classic()

ggplot() +
  geom_boxplot(data=taj.1a.ru.foc, aes(x = 1, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.1a.ru.back, aes(x = 2, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.10kb.ru, aes(x = 3, y = TajimaD),outlier.shape = NA) +
  geom_jitter(data=taj.1a.ru.foc, aes(x = 1, y = TajimaD),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=taj.1a.ru.rand, aes(x = 2, y = TajimaD),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=taj.all.ru.rand, aes(x = 3, y = TajimaD),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  ylim(-2.5,1) +
  theme_classic()

ggplot() +
  geom_boxplot(data=xpehh.1a.ruty.foc, aes(x = 1, y = XPEHH_rustica_tytleri),outlier.shape = NA) +
  geom_boxplot(data=xpehh.1a.ruty, aes(x = 2, y = XPEHH_rustica_tytleri),outlier.shape = NA) +
  geom_boxplot(data=xpehh.all.ruty, aes(x = 3, y = XPEHH_rustica_tytleri),outlier.shape = NA) +
  geom_jitter(data=xpehh.1a.ruty.foc, aes(x = 1, y = XPEHH_rustica_tytleri),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=xpehh.1a.ruty.rand, aes(x = 2, y = XPEHH_rustica_tytleri),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=xpehh.all.ruty.rand, aes(x = 3, y = XPEHH_rustica_tytleri),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  ylim(-5,10) +
  theme_classic()

## Also perform Mann-whitney tests between candidate region and backgrounds
wilcox.test(pbs.1a.ru.foc$pbs.10kb.ru,pbs.1a.ru.back$pbs.10kb.ru)
wilcox.test(pbs.1a.ru.foc$pbs.10kb.ru,pbs.10kb$pbs.10kb.ru)

wilcox.test(pi.1a.ru.foc$avg_pi,pi.1a.ru.back$avg_pi)
wilcox.test(pi.1a.ru.foc$avg_pi,pi.all.ru$avg_pi)

wilcox.test(taj.1a.ru.foc$TajimaD,taj.1a.ru.back$TajimaD)
wilcox.test(taj.1a.ru.foc$TajimaD,taj.10kb.ru$TajimaD)

wilcox.test(xpehh.1a.ruty.foc$XPEHH_rustica_tytleri,xpehh.1a.ruty$XPEHH_rustica_tytleri)
wilcox.test(xpehh.1a.ruty.foc$XPEHH_rustica_tytleri,xpehh.all.ruty$XPEHH_rustica_tytleri)

# Quick grab number of significant xp-EHH SNPs in KITLG association region (narrower than plotted region)
library(Rmisc)
xpehh.1a.ruty.foc <- xpehh.1a.ruty %>% filter(POSITION >= 4.415e+07 & POSITION <= 4.455e+07)
## Calculate bonferroni correction threshold
0.05/5928120
-log10(8.434377e-09)
count(xpehh.1a.ruty.foc$LOGPVALUE>=8.073947)

### Plot distributions of candidates and backgrounds (main figure)------------

# Set up backgrounds & candidate regions

## Fst
fst.all <- read.table('../pixy/pixy.all.order.fst.10kb.txt',header=T)
fst.all.ruty <- fst.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop1,'RU') & str_detect(pop2,'TY'))
fst.all.rugu <- fst.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop1,'RU') & str_detect(pop2,'GU'))
fst.all.guty <- fst.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop1,'GU') & str_detect(pop2,'TY'))

fst.all.ruty.1a.foc <- fst.all.ruty %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
fst.all.ruty.2.foc <- fst.all.ruty %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
fst.all.ruty.z.foc <- fst.all.ruty %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                         window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                         window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                         window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                         window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                         window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
fst.all.ruty.foc <- rbind(fst.all.ruty.1a.foc, fst.all.ruty.2.foc, fst.all.ruty.z.foc)

fst.all.rugu.1a.foc <- fst.all.rugu %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
fst.all.rugu.2.foc <- fst.all.rugu %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
fst.all.rugu.z.foc <- fst.all.rugu %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                         window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                         window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                         window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                         window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                         window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
fst.all.rugu.foc <- rbind(fst.all.rugu.1a.foc, fst.all.rugu.2.foc, fst.all.rugu.z.foc)

fst.all.guty.1a.foc <- fst.all.guty %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
fst.all.guty.2.foc <- fst.all.guty %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
fst.all.guty.z.foc <- fst.all.guty %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                         window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                         window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                         window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                         window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                         window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
fst.all.guty.foc <- rbind(fst.all.guty.1a.foc, fst.all.guty.2.foc, fst.all.guty.z.foc)


fst.all.ruty.foc <- fst.all.ruty.foc[which(fst.all.ruty.foc$avg_wc_fst!='NA'),]
fst.all.rugu.foc <- fst.all.rugu.foc[which(fst.all.rugu.foc$avg_wc_fst!='NA'),]
fst.all.guty.foc <- fst.all.guty.foc[which(fst.all.guty.foc$avg_wc_fst!='NA'),]

mean(fst.all.ruty.foc$avg_wc_fst,na.rm=T)
mean(fst.all.rugu.foc$avg_wc_fst,na.rm=T)
mean(fst.all.guty.foc$avg_wc_fst,na.rm=T)
std <- function(x) sd(x)/sqrt(length(x))

std(fst.all.ruty.foc$avg_wc_fst)
std(fst.all.rugu.foc$avg_wc_fst)
std(fst.all.guty.foc$avg_wc_fst)



## PBS
pbs.all.rand <- pbs.10kb[sample(nrow(pbs.10kb), 1200), ]
pbs.1a.foc <- pbs.10kb %>% filter(str_detect(fst.10kb.ruty.chromosome, 'NC_053453.1') & (fst.10kb.ruty.window_pos_1 >= 4.375e+07 & fst.10kb.ruty.window_pos_2 <= 4.5e+07 | fst.10kb.ruty.window_pos_1 >= 4.62e+07 & fst.10kb.ruty.window_pos_2 <= 4.69e+07))
pbs.2.foc <- pbs.10kb %>% filter(str_detect(fst.10kb.ruty.chromosome, 'NC_053450.1') & (fst.10kb.ruty.window_pos_1 >= 9.78e+07 & fst.10kb.ruty.window_pos_2 <= 9.855e+07 | fst.10kb.ruty.window_pos_1 >= 1.e+08 & fst.10kb.ruty.window_pos_2 <= 1.015e+08))
pbs.z.foc <- pbs.10kb %>% filter(str_detect(fst.10kb.ruty.chromosome, 'NC_053488.1') & (fst.10kb.ruty.window_pos_1 >= 1.85e+07 & fst.10kb.ruty.window_pos_2 <= 1.953e+07 | 
                                                                                          fst.10kb.ruty.window_pos_1 >= 3.275e+07 & fst.10kb.ruty.window_pos_2 <= 3.375e+07 | 
                                                                                          fst.10kb.ruty.window_pos_1 >= 3.5e+07 & fst.10kb.ruty.window_pos_2 <= 3.69e+07 | 
                                                                                          fst.10kb.ruty.window_pos_1 >= 3.5e+07 & fst.10kb.ruty.window_pos_2 <= 3.69e+07 | 
                                                                                          fst.10kb.ruty.window_pos_1 >= 4.36e+07 & fst.10kb.ruty.window_pos_2 <= 4.45e+07 | 
                                                                                          fst.10kb.ruty.window_pos_1 >= 4.46e+07 & fst.10kb.ruty.window_pos_2 <= 4.83e+07))

pbs.foc <- rbind(pbs.1a.foc,pbs.2.foc,pbs.z.foc)

## pi
pi.all <- read.table('../pixy/pixy.all.order.pi.10kb.txt',header=T)
pi.all.ru <- pi.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop,'RU'))
pi.all.ty <- pi.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop,'TY'))
pi.all.gu <- pi.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop,'GU'))
pi.all.ru.rand <- pi.all.ru[sample(nrow(pi.all.ru), 1200),]
pi.all.ty.rand <- pi.all.ty[sample(nrow(pi.all.ty), 1200),]
pi.all.gu.rand <- pi.all.gu[sample(nrow(pi.all.gu), 1200),]

pi.all.ru.1a.foc <- pi.all.ru %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
pi.all.ru.2.foc <- pi.all.ru %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
pi.all.ru.z.foc <- pi.all.ru %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                          window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                          window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                          window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                          window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                          window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
pi.all.ru.foc <- rbind(pi.all.ru.1a.foc, pi.all.ru.2.foc, pi.all.ru.z.foc)

pi.all.ty.1a.foc <- pi.all.ty %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
pi.all.ty.2.foc <- pi.all.ty %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
pi.all.ty.z.foc <- pi.all.ty %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                   window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                   window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
pi.all.ty.foc <- rbind(pi.all.ty.1a.foc, pi.all.ty.2.foc, pi.all.ty.z.foc)

pi.all.gu.1a.foc <- pi.all.gu %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
pi.all.gu.2.foc <- pi.all.gu %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
pi.all.gu.z.foc <- pi.all.gu %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                   window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                   window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
pi.all.gu.foc <- rbind(pi.all.gu.1a.foc, pi.all.gu.2.foc, pi.all.gu.z.foc)

## dxy
dxy.all <- read.table('../pixy/pixy.all.order.dxy.10kb.txt',header=T)
dxy.all.ruty <- dxy.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop1,'RU') & str_detect(pop2,'TY'))
dxy.all.rugu <- dxy.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop1,'RU') & str_detect(pop2,'GU'))
dxy.all.guty <- dxy.all %>% filter(!str_detect(chromosome,'NC_053487.1') & str_detect(pop1,'GU') & str_detect(pop2,'TY'))
dxy.all.ruty.rand <- dxy.all.ruty[sample(nrow(dxy.all.ruty), 1200),]
dxy.all.rugu.rand <- dxy.all.rugu[sample(nrow(dxy.all.rugu), 1200),]
dxy.all.guty.rand <- dxy.all.guty[sample(nrow(dxy.all.guty), 1200),]

dxy.all.ruty.1a.foc <- dxy.all.ruty %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
dxy.all.ruty.2.foc <- dxy.all.ruty %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
dxy.all.ruty.z.foc <- dxy.all.ruty %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                   window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                   window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
dxy.all.ruty.foc <- rbind(dxy.all.ruty.1a.foc, dxy.all.ruty.2.foc, dxy.all.ruty.z.foc)

dxy.all.rugu.1a.foc <- dxy.all.rugu %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
dxy.all.rugu.2.foc <- dxy.all.rugu %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
dxy.all.rugu.z.foc <- dxy.all.rugu %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                   window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                   window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
dxy.all.rugu.foc <- rbind(dxy.all.rugu.1a.foc, dxy.all.rugu.2.foc, dxy.all.rugu.z.foc)

dxy.all.guty.1a.foc <- dxy.all.guty %>% filter(str_detect(chromosome, 'NC_053453.1') & (window_pos_1 >= 4.375e+07 & window_pos_2 <= 4.5e+07 | window_pos_1 >= 4.62e+07 & window_pos_2 <= 4.69e+07))
dxy.all.guty.2.foc <- dxy.all.guty %>% filter(str_detect(chromosome, 'NC_053450.1') & (window_pos_1 >= 9.78e+07 & window_pos_2 <= 9.855e+07 | window_pos_1 >= 1.e+08 & window_pos_2 <= 1.015e+08))
dxy.all.guty.z.foc <- dxy.all.guty %>% filter(str_detect(chromosome, 'NC_053488.1') & (window_pos_1 >= 1.85e+07 & window_pos_2 <= 1.953e+07 | 
                                                                                   window_pos_1 >= 3.275e+07 & window_pos_2 <= 3.375e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 3.5e+07 & window_pos_2 <= 3.69e+07 | 
                                                                                   window_pos_1 >= 4.36e+07 & window_pos_2 <= 4.45e+07 | 
                                                                                   window_pos_1 >= 4.46e+07 & window_pos_2 <= 4.83e+07))
dxy.all.guty.foc <- rbind(dxy.all.guty.1a.foc, dxy.all.guty.2.foc, dxy.all.guty.z.foc)


## Tajima's D
taj.10kb.ru.rand <- taj.10kb.ru[sample(nrow(taj.10kb.ru), 1200),]
taj.10kb.ty.rand <- taj.10kb.ty[sample(nrow(taj.10kb.ty), 1200),]
taj.10kb.gu.rand <- taj.10kb.gu[sample(nrow(taj.10kb.gu), 1200),]

taj.10kb.ru.1a.foc <- taj.10kb.ru %>% filter(str_detect(CHROM, 'NC_053453.1') & (BIN_START >= 4.375e+07 & BIN_END <= 4.5e+07 | BIN_START >= 4.62e+07 & BIN_END <= 4.69e+07))
taj.10kb.ru.2.foc <- taj.10kb.ru %>% filter(str_detect(CHROM, 'NC_053450.1') & (BIN_START >= 9.78e+07 & BIN_END <= 9.855e+07 | BIN_START >= 1.e+08 & BIN_END <= 1.015e+08))
taj.10kb.ru.z.foc <- taj.10kb.ru %>% filter(str_detect(CHROM, 'NC_053488.1') & (BIN_START >= 1.85e+07 & BIN_END <= 1.953e+07 | 
                                                                                         BIN_START >= 3.275e+07 & BIN_END <= 3.375e+07 | 
                                                                                         BIN_START >= 3.5e+07 & BIN_END <= 3.69e+07 | 
                                                                                         BIN_START >= 3.5e+07 & BIN_END <= 3.69e+07 | 
                                                                                         BIN_START >= 4.36e+07 & BIN_END <= 4.45e+07 | 
                                                                                         BIN_START >= 4.46e+07 & BIN_END <= 4.83e+07))
taj.10kb.ru.foc <- rbind(taj.10kb.ru.1a.foc, taj.10kb.ru.2.foc, taj.10kb.ru.z.foc)

taj.10kb.ty.1a.foc <- taj.10kb.ty %>% filter(str_detect(CHROM, 'NC_053453.1') & (BIN_START >= 4.375e+07 & BIN_END <= 4.5e+07 | BIN_START >= 4.62e+07 & BIN_END <= 4.69e+07))
taj.10kb.ty.2.foc <- taj.10kb.ty %>% filter(str_detect(CHROM, 'NC_053450.1') & (BIN_START >= 9.78e+07 & BIN_END <= 9.855e+07 | BIN_START >= 1.e+08 & BIN_END <= 1.015e+08))
taj.10kb.ty.z.foc <- taj.10kb.ty %>% filter(str_detect(CHROM, 'NC_053488.1') & (BIN_START >= 1.85e+07 & BIN_END <= 1.953e+07 | 
                                                                                  BIN_START >= 3.275e+07 & BIN_END <= 3.375e+07 | 
                                                                                  BIN_START >= 3.5e+07 & BIN_END <= 3.69e+07 | 
                                                                                  BIN_START >= 3.5e+07 & BIN_END <= 3.69e+07 | 
                                                                                  BIN_START >= 4.36e+07 & BIN_END <= 4.45e+07 | 
                                                                                  BIN_START >= 4.46e+07 & BIN_END <= 4.83e+07))
taj.10kb.ty.foc <- rbind(taj.10kb.ty.1a.foc, taj.10kb.ty.2.foc, taj.10kb.ty.z.foc)

taj.10kb.gu.1a.foc <- taj.10kb.gu %>% filter(str_detect(CHROM, 'NC_053453.1') & (BIN_START >= 4.375e+07 & BIN_END <= 4.5e+07 | BIN_START >= 4.62e+07 & BIN_END <= 4.69e+07))
taj.10kb.gu.2.foc <- taj.10kb.gu %>% filter(str_detect(CHROM, 'NC_053450.1') & (BIN_START >= 9.78e+07 & BIN_END <= 9.855e+07 | BIN_START >= 1.e+08 & BIN_END <= 1.015e+08))
taj.10kb.gu.z.foc <- taj.10kb.gu %>% filter(str_detect(CHROM, 'NC_053488.1') & (BIN_START >= 1.85e+07 & BIN_END <= 1.953e+07 | 
                                                                                  BIN_START >= 3.275e+07 & BIN_END <= 3.375e+07 | 
                                                                                  BIN_START >= 3.5e+07 & BIN_END <= 3.69e+07 | 
                                                                                  BIN_START >= 3.5e+07 & BIN_END <= 3.69e+07 | 
                                                                                  BIN_START >= 4.36e+07 & BIN_END <= 4.45e+07 | 
                                                                                  BIN_START >= 4.46e+07 & BIN_END <= 4.83e+07))
taj.10kb.gu.foc <- rbind(taj.10kb.gu.1a.foc, taj.10kb.gu.2.foc, taj.10kb.gu.z.foc)


## iHS
ihs.ru.all.rand <- ihs.all.ru[sample(nrow(ihs.all.ru), 1200), ]
ihs.ty.all.rand <- ihs.all.ty[sample(nrow(ihs.all.ty), 1200), ]
ihs.gu.all.rand <- ihs.all.gu[sample(nrow(ihs.all.gu), 1200), ]

ihs.ru.1a.foc <- ihs.all.ru %>% filter(str_detect(CHR, 'NC_053453.1') & (POSITION >= 4.375e+07 & POSITION <= 4.5e+07 | POSITION >= 4.62e+07 & POSITION <= 4.69e+07))
ihs.ru.2.foc <- ihs.all.ru %>% filter(str_detect(CHR, 'NC_053450.1') & (POSITION >= 9.78e+07 & POSITION <= 9.855e+07 | POSITION >= 1.e+08 & POSITION <= 1.015e+08))
ihs.ru.z.foc <- ihs.all.ru %>% filter(str_detect(CHR, 'NC_053488.1') & (POSITION >= 1.85e+07 & POSITION <= 1.953e+07 | 
                                                                                          POSITION >= 3.275e+07 & POSITION <= 3.375e+07 | 
                                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                                          POSITION >= 4.36e+07 & POSITION <= 4.45e+07 | 
                                                                                          POSITION >= 4.46e+07 & POSITION <= 4.83e+07))
ihs.ru.foc <- rbind(ihs.ru.1a.foc, ihs.ru.2.foc, ihs.ru.z.foc)
ihs.ru.foc.rand <- ihs.ru.foc[sample(nrow(ihs.ru.foc), 1200),]

ihs.ty.1a.foc <- ihs.all.ty %>% filter(str_detect(CHR, 'NC_053453.1') & (POSITION >= 4.375e+07 & POSITION <= 4.5e+07 | POSITION >= 4.62e+07 & POSITION <= 4.69e+07))
ihs.ty.2.foc <- ihs.all.ty %>% filter(str_detect(CHR, 'NC_053450.1') & (POSITION >= 9.78e+07 & POSITION <= 9.855e+07 | POSITION >= 1.e+08 & POSITION <= 1.015e+08))
ihs.ty.z.foc <- ihs.all.ty %>% filter(str_detect(CHR, 'NC_053488.1') & (POSITION >= 1.85e+07 & POSITION <= 1.953e+07 | 
                                                                          POSITION >= 3.275e+07 & POSITION <= 3.375e+07 | 
                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                          POSITION >= 4.36e+07 & POSITION <= 4.45e+07 | 
                                                                          POSITION >= 4.46e+07 & POSITION <= 4.83e+07))
ihs.ty.foc <- rbind(ihs.ty.1a.foc, ihs.ty.2.foc, ihs.ty.z.foc)
ihs.ty.foc.rand <- ihs.ty.foc[sample(nrow(ihs.ty.foc), 1200),]

ihs.gu.1a.foc <- ihs.all.gu %>% filter(str_detect(CHR, 'NC_053453.1') & (POSITION >= 4.375e+07 & POSITION <= 4.5e+07 | POSITION >= 4.62e+07 & POSITION <= 4.69e+07))
ihs.gu.2.foc <- ihs.all.gu %>% filter(str_detect(CHR, 'NC_053450.1') & (POSITION >= 9.78e+07 & POSITION <= 9.855e+07 | POSITION >= 1.e+08 & POSITION <= 1.015e+08))
ihs.gu.z.foc <- ihs.all.gu %>% filter(str_detect(CHR, 'NC_053488.1') & (POSITION >= 1.85e+07 & POSITION <= 1.953e+07 | 
                                                                          POSITION >= 3.275e+07 & POSITION <= 3.375e+07 | 
                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                          POSITION >= 4.36e+07 & POSITION <= 4.45e+07 | 
                                                                          POSITION >= 4.46e+07 & POSITION <= 4.83e+07))
ihs.gu.foc <- rbind(ihs.gu.1a.foc, ihs.gu.2.foc, ihs.gu.z.foc)
ihs.gu.foc.rand <- ihs.gu.foc[sample(nrow(ihs.gu.foc), 1200),]


## xp-EHH
xpehh.ruty.all.rand <- xpehh.all.ruty[sample(nrow(xpehh.all.ruty), 1200), ]
xpehh.rugu.all.rand <- xpehh.all.rugu[sample(nrow(xpehh.all.rugu), 1200), ]
xpehh.guty.all.rand <- xpehh.all.guty[sample(nrow(xpehh.all.guty), 1200), ]

xpehh.ruty.1a.foc <- xpehh.all.ruty %>% filter(str_detect(CHR, 'NC_053453.1') & (POSITION >= 4.375e+07 & POSITION <= 4.5e+07 | POSITION >= 4.62e+07 & POSITION <= 4.69e+07))
xpehh.ruty.2.foc <- xpehh.all.ruty %>% filter(str_detect(CHR, 'NC_053450.1') & (POSITION >= 9.78e+07 & POSITION <= 9.855e+07 | POSITION >= 1.e+08 & POSITION <= 1.015e+08))
xpehh.ruty.z.foc <- xpehh.all.ruty %>% filter(str_detect(CHR, 'NC_053488.1') & (POSITION >= 1.85e+07 & POSITION <= 1.953e+07 | 
                                                                          POSITION >= 3.275e+07 & POSITION <= 3.375e+07 | 
                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                          POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                          POSITION >= 4.36e+07 & POSITION <= 4.45e+07 | 
                                                                          POSITION >= 4.46e+07 & POSITION <= 4.83e+07))
xpehh.ruty.foc <- rbind(xpehh.ruty.1a.foc, xpehh.ruty.2.foc, xpehh.ruty.z.foc)
xpehh.ruty.foc.rand <- xpehh.ruty.foc[sample(nrow(xpehh.ruty.foc), 1200),]

xpehh.rugu.1a.foc <- xpehh.all.rugu %>% filter(str_detect(CHR, 'NC_053453.1') & (POSITION >= 4.375e+07 & POSITION <= 4.5e+07 | POSITION >= 4.62e+07 & POSITION <= 4.69e+07))
xpehh.rugu.2.foc <- xpehh.all.rugu %>% filter(str_detect(CHR, 'NC_053450.1') & (POSITION >= 9.78e+07 & POSITION <= 9.855e+07 | POSITION >= 1.e+08 & POSITION <= 1.015e+08))
xpehh.rugu.z.foc <- xpehh.all.rugu %>% filter(str_detect(CHR, 'NC_053488.1') & (POSITION >= 1.85e+07 & POSITION <= 1.953e+07 | 
                                                                                  POSITION >= 3.275e+07 & POSITION <= 3.375e+07 | 
                                                                                  POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                                  POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                                  POSITION >= 4.36e+07 & POSITION <= 4.45e+07 | 
                                                                                  POSITION >= 4.46e+07 & POSITION <= 4.83e+07))
xpehh.rugu.foc <- rbind(xpehh.rugu.1a.foc, xpehh.rugu.2.foc, xpehh.rugu.z.foc)
xpehh.rugu.foc.rand <- xpehh.rugu.foc[sample(nrow(xpehh.rugu.foc), 1200),]

xpehh.guty.1a.foc <- xpehh.all.guty %>% filter(str_detect(CHR, 'NC_053453.1') & (POSITION >= 4.375e+07 & POSITION <= 4.5e+07 | POSITION >= 4.62e+07 & POSITION <= 4.69e+07))
xpehh.guty.2.foc <- xpehh.all.guty %>% filter(str_detect(CHR, 'NC_053450.1') & (POSITION >= 9.78e+07 & POSITION <= 9.855e+07 | POSITION >= 1.e+08 & POSITION <= 1.015e+08))
xpehh.guty.z.foc <- xpehh.all.guty %>% filter(str_detect(CHR, 'NC_053488.1') & (POSITION >= 1.85e+07 & POSITION <= 1.953e+07 | 
                                                                                  POSITION >= 3.275e+07 & POSITION <= 3.375e+07 | 
                                                                                  POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                                  POSITION >= 3.5e+07 & POSITION <= 3.69e+07 | 
                                                                                  POSITION >= 4.36e+07 & POSITION <= 4.45e+07 | 
                                                                                  POSITION >= 4.46e+07 & POSITION <= 4.83e+07))
xpehh.guty.foc <- rbind(xpehh.guty.1a.foc, xpehh.guty.2.foc, xpehh.guty.z.foc)
xpehh.guty.foc.rand <- xpehh.guty.foc[sample(nrow(xpehh.guty.foc), 1200),]

## Perform statistical tests comparing background and candidate distributions
wilcox.test(pbs.10kb$pbs.10kb.ru,pbs.foc$pbs.10kb.ru)
wilcox.test(pbs.10kb$pbs.10kb.ty,pbs.foc$pbs.10kb.ty)
wilcox.test(pbs.10kb$pbs.10kb.gu,pbs.foc$pbs.10kb.gu)

wilcox.test(pi.all.ru$avg_pi,pi.all.ru.foc$avg_pi)
wilcox.test(pi.all.ty$avg_pi,pi.all.ty.foc$avg_pi)
wilcox.test(pi.all.gu$avg_pi,pi.all.gu.foc$avg_pi)

wilcox.test(abs(ihs.all.ru$IHS),abs(ihs.ru.foc$IHS))
wilcox.test(abs(ihs.all.ty$IHS),abs(ihs.ty.foc$IHS))
wilcox.test(abs(ihs.all.gu$IHS),abs(ihs.gu.foc$IHS))

ks.test(xpehh.all.ruty$XPEHH_rustica_tytleri,xpehh.ruty.foc$XPEHH_rustica_tytleri)
ks.test(xpehh.all.rugu$XPEHH_rustica_gutturalis,xpehh.rugu.foc$XPEHH_rustica_gutturalis)
ks.test(xpehh.all.guty$XPEHH_tytleri_gutturalis,xpehh.guty.foc$XPEHH_tytleri_gutturalis)

## Plots 
ggplot() +
  geom_boxplot(data=pbs.10kb, aes(x = 1, y = pbs.10kb.ru),outlier.shape = NA) +
  geom_boxplot(data=pbs.foc, aes(x = 2, y = pbs.10kb.ru),outlier.shape = NA) +
  geom_boxplot(data=pbs.10kb, aes(x = 3, y = pbs.10kb.ty),outlier.shape = NA) +
  geom_boxplot(data=pbs.foc, aes(x = 4, y = pbs.10kb.ty),outlier.shape = NA) +
  geom_boxplot(data=pbs.10kb, aes(x = 5, y = pbs.10kb.gu),outlier.shape = NA) +
  geom_boxplot(data=pbs.foc, aes(x = 6, y = pbs.10kb.gu),outlier.shape = NA) +
  geom_jitter(data=pbs.all.rand, aes(x = 1, y = pbs.10kb.ru),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=pbs.foc, aes(x = 2, y = pbs.10kb.ru),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=pbs.all.rand, aes(x = 3, y = pbs.10kb.ty),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=pbs.foc, aes(x = 4, y = pbs.10kb.ty),width = 0.075,colour = alpha('#EBB320',0.25),cex=0.5) +
  geom_jitter(data=pbs.all.rand, aes(x = 5, y = pbs.10kb.gu),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=pbs.foc, aes(x = 6, y = pbs.10kb.gu),width = 0.075,colour = alpha('#6BA5CC',0.25),cex=0.5) +
  ylab('PBS')+
  ylim(0,1) +
  theme_classic()

ggplot() +
  geom_boxplot(data=pi.all.ru, aes(x = 1, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.all.ru.foc, aes(x = 2, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.all.ty, aes(x = 3, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.all.ty.foc, aes(x = 4, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.all.gu, aes(x = 5, y = avg_pi),outlier.shape = NA) +
  geom_boxplot(data=pi.all.gu.foc, aes(x = 6, y = avg_pi),outlier.shape = NA) +
  geom_jitter(data=pi.all.ru.rand, aes(x = 1, y = avg_pi),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=pi.all.ru.foc, aes(x = 2, y = avg_pi),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=pi.all.ty.rand, aes(x = 3, y = avg_pi),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=pi.all.ty.foc, aes(x = 4, y = avg_pi),width = 0.075,colour = alpha('#EBB320',0.25),cex=0.5) +
  geom_jitter(data=pi.all.gu.rand, aes(x = 5, y = avg_pi),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=pi.all.gu.foc, aes(x = 6, y = avg_pi),width = 0.075,colour = alpha('#6BA5CC',0.25),cex=0.5) +
  ylab('π')+
  ylim(0,0.01) +
  theme_classic()

ggplot() +
  geom_boxplot(data=taj.10kb.ru, aes(x = 1, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.10kb.ru.foc, aes(x = 2, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.10kb.ty, aes(x = 3, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.10kb.ty.foc, aes(x = 4, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.10kb.gu, aes(x = 5, y = TajimaD),outlier.shape = NA) +
  geom_boxplot(data=taj.10kb.gu.foc, aes(x = 6, y = TajimaD),outlier.shape = NA) +
  geom_jitter(data=taj.10kb.ru.rand, aes(x = 1, y = TajimaD),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=taj.10kb.ru.foc, aes(x = 2, y = TajimaD),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=taj.10kb.ty.rand, aes(x = 3, y = TajimaD),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=taj.10kb.ty.foc, aes(x = 4, y = TajimaD),width = 0.075,colour = alpha('#EBB320',0.25),cex=0.5) +
  geom_jitter(data=taj.10kb.gu.rand, aes(x = 5, y = TajimaD),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=taj.10kb.gu.foc, aes(x = 6, y = TajimaD),width = 0.075,colour = alpha('#6BA5CC',0.25),cex=0.5) +
  ylab('π')+
  ylim(-3,3) +
  theme_classic()

ggplot() +
  geom_boxplot(data=dxy.all.ruty, aes(x = 1, y = avg_dxy),outlier.shape = NA) +
  geom_boxplot(data=dxy.all.ruty.foc, aes(x = 2, y = avg_dxy),outlier.shape = NA) +
  geom_boxplot(data=dxy.all.rugu, aes(x = 3, y = avg_dxy),outlier.shape = NA) +
  geom_boxplot(data=dxy.all.rugu.foc, aes(x = 4, y = avg_dxy),outlier.shape = NA) +
  geom_boxplot(data=dxy.all.guty, aes(x = 5, y = avg_dxy),outlier.shape = NA) +
  geom_boxplot(data=dxy.all.guty.foc, aes(x = 6, y = avg_dxy),outlier.shape = NA) +
  geom_jitter(data=dxy.all.ruty.rand, aes(x = 1, y = avg_dxy),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=dxy.all.ruty.foc, aes(x = 2, y = avg_dxy),width = 0.075,colour = alpha('#E28026',0.25),cex=0.5) +
  geom_jitter(data=dxy.all.rugu.rand, aes(x = 3, y = avg_dxy),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=dxy.all.rugu.foc, aes(x = 4, y = avg_dxy),width = 0.075,colour = alpha('#8362AA',0.25),cex=0.5) +
  geom_jitter(data=dxy.all.guty.rand, aes(x = 5, y = avg_dxy),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=dxy.all.guty.foc, aes(x = 6, y = avg_dxy),width = 0.075,colour = alpha('#70C6A8',0.25),cex=0.5) +
  ylab('dxy')+
  ylim(0,0.01) +
  theme_classic()

ggplot() +
  geom_boxplot(data=ihs.all.ru, aes(x = 1, y = abs(IHS)),outlier.shape = NA) +
  geom_boxplot(data=ihs.ru.foc, aes(x = 2, y = abs(IHS)),outlier.shape = NA) +
  geom_boxplot(data=ihs.all.ty, aes(x = 3, y = abs(IHS)),outlier.shape = NA) +
  geom_boxplot(data=ihs.ty.foc, aes(x = 4, y = abs(IHS)),outlier.shape = NA) +
  geom_boxplot(data=ihs.all.gu, aes(x = 5, y = abs(IHS)),outlier.shape = NA) +
  geom_boxplot(data=ihs.gu.foc, aes(x = 6, y = abs(IHS)),outlier.shape = NA) +
  geom_jitter(data=ihs.ru.all.rand, aes(x = 1, y = abs(IHS)),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=ihs.ru.foc.rand, aes(x = 2, y = abs(IHS)),width = 0.075,colour = alpha('#B03160',0.25),cex=0.5) +
  geom_jitter(data=ihs.ty.all.rand, aes(x = 3, y = abs(IHS)),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=ihs.ty.foc.rand, aes(x = 4, y = abs(IHS)),width = 0.075,colour = alpha('#EBB320',0.25),cex=0.5) +
  geom_jitter(data=ihs.gu.all.rand, aes(x = 5, y = abs(IHS)),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=ihs.gu.foc.rand, aes(x = 6, y = abs(IHS)),width = 0.075,colour = alpha('#6BA5CC',0.25),cex=0.5) +
  ylab('|iHS|')+
  ylim(0,4) +
  theme_classic()

ggplot() +
  geom_boxplot(data=xpehh.all.ruty, aes(x = 1, y = XPEHH_rustica_tytleri),outlier.shape = NA) +
  geom_boxplot(data=xpehh.ruty.foc, aes(x = 2, y = XPEHH_rustica_tytleri),outlier.shape = NA) +
  geom_boxplot(data=xpehh.all.rugu, aes(x = 3, y = XPEHH_rustica_gutturalis),outlier.shape = NA) +
  geom_boxplot(data=xpehh.rugu.foc, aes(x = 4, y = XPEHH_rustica_gutturalis),outlier.shape = NA) +
  geom_boxplot(data=xpehh.all.guty, aes(x = 5, y = XPEHH_tytleri_gutturalis),outlier.shape = NA) +
  geom_boxplot(data=xpehh.guty.foc, aes(x = 6, y = XPEHH_tytleri_gutturalis),outlier.shape = NA) +
  geom_jitter(data=xpehh.ruty.all.rand, aes(x = 1, y = XPEHH_rustica_tytleri),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=xpehh.ruty.foc.rand, aes(x = 2, y = XPEHH_rustica_tytleri),width = 0.075,colour = alpha('#E28026',0.25),cex=0.5) +
  geom_jitter(data=xpehh.rugu.all.rand, aes(x = 3, y = XPEHH_rustica_gutturalis),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=xpehh.rugu.foc.rand, aes(x = 4, y = XPEHH_rustica_gutturalis),width = 0.075,colour = alpha('#8362AA',0.25),cex=0.5) +
  geom_jitter(data=xpehh.guty.all.rand, aes(x = 5, y = XPEHH_tytleri_gutturalis),width = 0.075,colour = alpha('grey',0.25),cex=0.5) +
  geom_jitter(data=xpehh.guty.foc.rand, aes(x = 6, y = XPEHH_tytleri_gutturalis),width = 0.075,colour = alpha('#70C6A8',0.25),cex=0.5) +
  ylab('xp-EHH')+
  ylim(-5,5) +
  theme_classic()

### Plot candidate regions----------------------------------------------------

# Output at 6 x 5.5 (may change this; this is from multi-panel GWA + Fst highlights)
# Export at 12 x 3.5

## KITLG region
x <- c(4.375e+07,4.5e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053453.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053453.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.1a.ruty.foc <- xpehh.1a.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.1a.rugu.foc <- xpehh.1a.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.1a.guty.foc <- xpehh.1a.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.1a.ruty.foc$color <- 'lightgrey'
xpehh.1a.ruty.foc$color[xpehh.1a.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.1a.ruty.foc$color[xpehh.1a.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.1a.rugu.foc$color <- 'lightgrey'
xpehh.1a.rugu.foc$color[xpehh.1a.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.1a.rugu.foc$color[xpehh.1a.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.1a.guty.foc$color <- 'lightgrey'
xpehh.1a.guty.foc$color[xpehh.1a.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.1a.guty.foc$color[xpehh.1a.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.50
y2 <- y-0.1
plot(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.ru,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.6))
lines(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.1a.ru$window_pos_1+(pi.50kb.1a.ru$window_pos_2 - pi.50kb.1a.ru$window_pos_1)/2,pi.50kb.1a.ru$avg_pi,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.0025))
lines(pi.50kb.1a.ty$window_pos_1+(pi.50kb.1a.ty$window_pos_2 - pi.50kb.1a.ty$window_pos_1)/2,pi.50kb.1a.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.1a.gu$window_pos_1+(pi.50kb.1a.gu$window_pos_2 - pi.50kb.1a.gu$window_pos_1)/2,pi.50kb.1a.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.1a.ru$BIN_START+(taj.50kb.1a.ru$BIN_END - taj.50kb.1a.ru$BIN_START)/2,taj.50kb.1a.ru$TajimaD,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.5,1))
lines(taj.50kb.1a.ty$BIN_START+(taj.50kb.1a.ty$BIN_END - taj.50kb.1a.ty$BIN_START)/2,taj.50kb.1a.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.1a.gu$BIN_START+(taj.50kb.1a.gu$BIN_END - taj.50kb.1a.gu$BIN_START)/2,taj.50kb.1a.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.1a.ruty.foc$POSITION,xpehh.1a.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.1a.ruty.foc$color,xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.1a.rugu.foc$POSITION,xpehh.1a.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.1a.rugu.foc$color,xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.1a.guty.foc$POSITION,xpehh.1a.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.1a.guty.foc$color,xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## PLXNC1 region
x <- c(4.62e+07,4.69e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053453.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053453.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053453.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.1a.ruty.foc <- xpehh.1a.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.1a.rugu.foc <- xpehh.1a.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.1a.guty.foc <- xpehh.1a.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.1a.ruty.foc$color <- 'lightgrey'
xpehh.1a.ruty.foc$color[xpehh.1a.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.1a.ruty.foc$color[xpehh.1a.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.1a.rugu.foc$color <- 'lightgrey'
xpehh.1a.rugu.foc$color[xpehh.1a.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.1a.rugu.foc$color[xpehh.1a.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.1a.guty.foc$color <- 'lightgrey'
xpehh.1a.guty.foc$color[xpehh.1a.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.1a.guty.foc$color[xpehh.1a.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.50
y2 <- y-0.1
plot(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.ru,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.6))
lines(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1+(pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_2-pbs.50kb.1a$fst.50kb.1a.ruty.window_pos_1),pbs.50kb.1a$pbs.50kb.1a.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.1a.ru$window_pos_1+(pi.50kb.1a.ru$window_pos_2 - pi.50kb.1a.ru$window_pos_1)/2,pi.50kb.1a.ru$avg_pi,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.007))
lines(pi.50kb.1a.ty$window_pos_1+(pi.50kb.1a.ty$window_pos_2 - pi.50kb.1a.ty$window_pos_1)/2,pi.50kb.1a.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.1a.gu$window_pos_1+(pi.50kb.1a.gu$window_pos_2 - pi.50kb.1a.gu$window_pos_1)/2,pi.50kb.1a.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.1a.ru$BIN_START+(taj.50kb.1a.ru$BIN_END - taj.50kb.1a.ru$BIN_START)/2,taj.50kb.1a.ru$TajimaD,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2,0))
lines(taj.50kb.1a.ty$BIN_START+(taj.50kb.1a.ty$BIN_END - taj.50kb.1a.ty$BIN_START)/2,taj.50kb.1a.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.1a.gu$BIN_START+(taj.50kb.1a.gu$BIN_END - taj.50kb.1a.gu$BIN_START)/2,taj.50kb.1a.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.1a.ruty.foc$POSITION,xpehh.1a.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.1a.ruty.foc$color,xlim=x,ylim=c(-7,7),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.1a.rugu.foc$POSITION,xpehh.1a.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.1a.rugu.foc$color,xlim=x,ylim=c(-7,7),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.1a.guty.foc$POSITION,xpehh.1a.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.1a.guty.foc$color,xlim=x,ylim=c(-7,7),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## SPEF2/PRLR/SLC45A2 region
x <- c(1.85e+07,1.953e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.z.ruty.foc <- xpehh.z.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.rugu.foc <- xpehh.z.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.guty.foc <- xpehh.z.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.ruty.foc$color <- 'lightgrey'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.z.rugu.foc$color <- 'lightgrey'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.z.guty.foc$color <- 'lightgrey'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.6
y2 <- y-0.1
plot(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ru,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.7))
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.z.ru$window_pos_1+(pi.50kb.z.ru$window_pos_2 - pi.50kb.z.ru$window_pos_1)/2,pi.50kb.z.ru$avg_pi,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.007))
lines(pi.50kb.z.ty$window_pos_1+(pi.50kb.z.ty$window_pos_2 - pi.50kb.z.ty$window_pos_1)/2,pi.50kb.z.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.z.gu$window_pos_1+(pi.50kb.z.gu$window_pos_2 - pi.50kb.z.gu$window_pos_1)/2,pi.50kb.z.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.z.ru$BIN_START+(taj.50kb.z.ru$BIN_END - taj.50kb.z.ru$BIN_START)/2,taj.50kb.z.ru$TajimaD,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.25,0))
lines(taj.50kb.z.ty$BIN_START+(taj.50kb.z.ty$BIN_END - taj.50kb.z.ty$BIN_START)/2,taj.50kb.z.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.z.gu$BIN_START+(taj.50kb.z.gu$BIN_END - taj.50kb.z.gu$BIN_START)/2,taj.50kb.z.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.z.ruty.foc$POSITION,xpehh.z.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.z.ruty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.rugu.foc$POSITION,xpehh.z.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.z.rugu.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.guty.foc$POSITION,xpehh.z.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.z.guty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## BNC2 region
x <- c(3.275e+07,3.375e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.z.ruty.foc <- xpehh.z.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.rugu.foc <- xpehh.z.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.guty.foc <- xpehh.z.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.ruty.foc$color <- 'lightgrey'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.z.rugu.foc$color <- 'lightgrey'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.z.guty.foc$color <- 'lightgrey'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.8
y2 <- y-0.1
plot(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ru,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.9))
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.z.ru$window_pos_1+(pi.50kb.z.ru$window_pos_2 - pi.50kb.z.ru$window_pos_1)/2,pi.50kb.z.ru$avg_pi,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.003))
lines(pi.50kb.z.ty$window_pos_1+(pi.50kb.z.ty$window_pos_2 - pi.50kb.z.ty$window_pos_1)/2,pi.50kb.z.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.z.gu$window_pos_1+(pi.50kb.z.gu$window_pos_2 - pi.50kb.z.gu$window_pos_1)/2,pi.50kb.z.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.z.ru$BIN_START+(taj.50kb.z.ru$BIN_END - taj.50kb.z.ru$BIN_START)/2,taj.50kb.z.ru$TajimaD,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.25,0))
lines(taj.50kb.z.ty$BIN_START+(taj.50kb.z.ty$BIN_END - taj.50kb.z.ty$BIN_START)/2,taj.50kb.z.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.z.gu$BIN_START+(taj.50kb.z.gu$BIN_END - taj.50kb.z.gu$BIN_START)/2,taj.50kb.z.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.z.ruty.foc$POSITION,xpehh.z.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.z.ruty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.rugu.foc$POSITION,xpehh.z.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.z.rugu.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.guty.foc$POSITION,xpehh.z.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.z.guty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## GNAQ region
x <- c(3.5e+07,3.69e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.z.ruty.foc <- xpehh.z.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.rugu.foc <- xpehh.z.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.guty.foc <- xpehh.z.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.ruty.foc$color <- 'lightgrey'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.z.rugu.foc$color <- 'lightgrey'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.z.guty.foc$color <- 'lightgrey'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.9
y2 <- y-0.2
plot(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ru,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,1.1))
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.z.ru$window_pos_1+(pi.50kb.z.ru$window_pos_2 - pi.50kb.z.ru$window_pos_1)/2,pi.50kb.z.ru$avg_pi,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.003))
lines(pi.50kb.z.ty$window_pos_1+(pi.50kb.z.ty$window_pos_2 - pi.50kb.z.ty$window_pos_1)/2,pi.50kb.z.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.z.gu$window_pos_1+(pi.50kb.z.gu$window_pos_2 - pi.50kb.z.gu$window_pos_1)/2,pi.50kb.z.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.z.ru$BIN_START+(taj.50kb.z.ru$BIN_END - taj.50kb.z.ru$BIN_START)/2,taj.50kb.z.ru$TajimaD,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.25,0))
lines(taj.50kb.z.ty$BIN_START+(taj.50kb.z.ty$BIN_END - taj.50kb.z.ty$BIN_START)/2,taj.50kb.z.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.z.gu$BIN_START+(taj.50kb.z.gu$BIN_END - taj.50kb.z.gu$BIN_START)/2,taj.50kb.z.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.z.ruty.foc$POSITION,xpehh.z.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.z.ruty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.rugu.foc$POSITION,xpehh.z.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.z.rugu.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.guty.foc$POSITION,xpehh.z.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.z.guty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## ROR2 region
x <- c(4.36e+07,4.45e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.z.ruty.foc <- xpehh.z.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.rugu.foc <- xpehh.z.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.guty.foc <- xpehh.z.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.ruty.foc$color <- 'lightgrey'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.z.rugu.foc$color <- 'lightgrey'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.z.guty.foc$color <- 'lightgrey'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 1.3
y2 <- y-0.2
plot(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ru,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,1.5))
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.z.ru$window_pos_1+(pi.50kb.z.ru$window_pos_2 - pi.50kb.z.ru$window_pos_1)/2,pi.50kb.z.ru$avg_pi,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.002))
lines(pi.50kb.z.ty$window_pos_1+(pi.50kb.z.ty$window_pos_2 - pi.50kb.z.ty$window_pos_1)/2,pi.50kb.z.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.z.gu$window_pos_1+(pi.50kb.z.gu$window_pos_2 - pi.50kb.z.gu$window_pos_1)/2,pi.50kb.z.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.z.ru$BIN_START+(taj.50kb.z.ru$BIN_END - taj.50kb.z.ru$BIN_START)/2,taj.50kb.z.ru$TajimaD,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.5,0))
lines(taj.50kb.z.ty$BIN_START+(taj.50kb.z.ty$BIN_END - taj.50kb.z.ty$BIN_START)/2,taj.50kb.z.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.z.gu$BIN_START+(taj.50kb.z.gu$BIN_END - taj.50kb.z.gu$BIN_START)/2,taj.50kb.z.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.z.ruty.foc$POSITION,xpehh.z.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.z.ruty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.rugu.foc$POSITION,xpehh.z.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.z.rugu.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.guty.foc$POSITION,xpehh.z.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.z.guty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## APC/CAMK4 region
x <- c(4.45e+07,4.83e+07)
cand.vc.for.foc <- cand.vc.for %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.vc.rev.foc <- cand.vc.rev %>% filter(str_detect(chrom.gene, 'NC_053488.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053488.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.z.ruty.foc <- xpehh.z.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.rugu.foc <- xpehh.z.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.guty.foc <- xpehh.z.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.z.ruty.foc$color <- 'lightgrey'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.z.ruty.foc$color[xpehh.z.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.z.rugu.foc$color <- 'lightgrey'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.z.rugu.foc$color[xpehh.z.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.z.guty.foc$color <- 'lightgrey'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.z.guty.foc$color[xpehh.z.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 1.4
y2 <- y-0.2
plot(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ru,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,1.5))
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.z$fst.50kb.z.ruty.window_pos_1+(pbs.50kb.z$fst.50kb.z.ruty.window_pos_2-pbs.50kb.z$fst.50kb.z.ruty.window_pos_1),pbs.50kb.z$pbs.50kb.z.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.05,xright=exons.for.foc$end,ybottom=y-0.05,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.05,xright=exons.rev.foc$end,ybottom=y2-0.05,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.vc.for.foc){arrows(x0=cand.vc.for.foc$start,y0=y,x1=cand.vc.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.vc.for.foc$end,y,cand.vc.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.vc.rev.foc){arrows(x0=cand.vc.rev.foc$start,y0=y2,x1=cand.vc.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.vc.rev.foc$end,y2,cand.vc.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.z.ru$window_pos_1+(pi.50kb.z.ru$window_pos_2 - pi.50kb.z.ru$window_pos_1)/2,pi.50kb.z.ru$avg_pi,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.001))
lines(pi.50kb.z.ty$window_pos_1+(pi.50kb.z.ty$window_pos_2 - pi.50kb.z.ty$window_pos_1)/2,pi.50kb.z.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.z.gu$window_pos_1+(pi.50kb.z.gu$window_pos_2 - pi.50kb.z.gu$window_pos_1)/2,pi.50kb.z.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.z.ru$BIN_START+(taj.50kb.z.ru$BIN_END - taj.50kb.z.ru$BIN_START)/2,taj.50kb.z.ru$TajimaD,type='l',lwd=1.5,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.5,0))
lines(taj.50kb.z.ty$BIN_START+(taj.50kb.z.ty$BIN_END - taj.50kb.z.ty$BIN_START)/2,taj.50kb.z.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.z.gu$BIN_START+(taj.50kb.z.gu$BIN_END - taj.50kb.z.gu$BIN_START)/2,taj.50kb.z.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.z.ruty.foc$POSITION,xpehh.z.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.z.ruty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.rugu.foc$POSITION,xpehh.z.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.z.rugu.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.z.guty.foc$POSITION,xpehh.z.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.z.guty.foc$color,xlim=x,ylim=c(-6,6),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## ICE1/lncRNA region
x <- c(9.78e+07,9.855e+07)
cand.ts.for.foc <- cand.ts.for %>% filter(str_detect(chrom.gene, 'NC_053450.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.ts.rev.foc <- cand.ts.rev %>% filter(str_detect(chrom.gene, 'NC_053450.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.2.ruty.foc <- xpehh.2.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.2.rugu.foc <- xpehh.2.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.2.guty.foc <- xpehh.2.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.2.ruty.foc$color <- 'lightgrey'
xpehh.2.ruty.foc$color[xpehh.2.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.2.ruty.foc$color[xpehh.2.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.2.rugu.foc$color <- 'lightgrey'
xpehh.2.rugu.foc$color[xpehh.2.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.2.rugu.foc$color[xpehh.2.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.2.guty.foc$color <- 'lightgrey'
xpehh.2.guty.foc$color[xpehh.2.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.2.guty.foc$color[xpehh.2.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.50
y2 <- y-0.1
plot(pbs.50kb.2$fst.50kb.2.ruty.window_pos_1+(pbs.50kb.2$fst.50kb.2.ruty.window_pos_2-pbs.50kb.2$fst.50kb.2.ruty.window_pos_1),pbs.50kb.2$pbs.50kb.2.ru,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.6))
lines(pbs.50kb.2$fst.50kb.2.ruty.window_pos_1+(pbs.50kb.2$fst.50kb.2.ruty.window_pos_2-pbs.50kb.2$fst.50kb.2.ruty.window_pos_1),pbs.50kb.2$pbs.50kb.2.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.2$fst.50kb.2.ruty.window_pos_1+(pbs.50kb.2$fst.50kb.2.ruty.window_pos_2-pbs.50kb.2$fst.50kb.2.ruty.window_pos_1),pbs.50kb.2$pbs.50kb.2.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.ts.for.foc){arrows(x0=cand.ts.for.foc$start,y0=y,x1=cand.ts.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.ts.for.foc$end,y,cand.ts.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.ts.rev.foc){arrows(x0=cand.ts.rev.foc$start,y0=y2,x1=cand.ts.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.ts.rev.foc$end,y2,cand.ts.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.2.ru$window_pos_1+(pi.50kb.2.ru$window_pos_2 - pi.50kb.2.ru$window_pos_1)/2,pi.50kb.2.ru$avg_pi,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.0025))
lines(pi.50kb.2.ty$window_pos_1+(pi.50kb.2.ty$window_pos_2 - pi.50kb.2.ty$window_pos_1)/2,pi.50kb.2.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.2.gu$window_pos_1+(pi.50kb.2.gu$window_pos_2 - pi.50kb.2.gu$window_pos_1)/2,pi.50kb.2.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.2.ru$BIN_START+(taj.50kb.2.ru$BIN_END - taj.50kb.2.ru$BIN_START)/2,taj.50kb.2.ru$TajimaD,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.5,1))
lines(taj.50kb.2.ty$BIN_START+(taj.50kb.2.ty$BIN_END - taj.50kb.2.ty$BIN_START)/2,taj.50kb.2.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.2.gu$BIN_START+(taj.50kb.2.gu$BIN_END - taj.50kb.2.gu$BIN_START)/2,taj.50kb.2.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.2.ruty.foc$POSITION,xpehh.2.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.2.ruty.foc$color,xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.2.rugu.foc$POSITION,xpehh.2.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.2.rugu.foc$color,xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.2.guty.foc$POSITION,xpehh.2.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.2.guty.foc$color,xlim=x,ylim=c(-5,10),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')


## PDE1C region
x <- c(1.e+08,1.015e+08)
cand.ts.for.foc <- cand.ts.for %>% filter(str_detect(chrom.gene, 'NC_053450.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
cand.ts.rev.foc <- cand.ts.rev %>% filter(str_detect(chrom.gene, 'NC_053450.1')) %>% filter(start.gene >= x[1]-10000 && end.gene <= x[2]+10000)
genes.for.foc <- genes.for %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
genes.rev.foc <- genes.rev %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.for.foc <- exons.for %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
exons.rev.foc <- exons.rev %>% filter(str_detect(chrom, 'NC_053450.1')) %>% filter(start >= x[1]-10000 & end <= x[2]+10000)
xpehh.2.ruty.foc <- xpehh.2.ruty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.2.rugu.foc <- xpehh.2.rugu %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.2.guty.foc <- xpehh.2.guty %>% filter(POSITION >= x[1]-100000 & POSITION <= x[2]+100000)
xpehh.2.ruty.foc$color <- 'lightgrey'
xpehh.2.ruty.foc$color[xpehh.2.ruty.foc$XPEHH_rustica_tytleri < -2.497729] <- '#EBB320'
xpehh.2.ruty.foc$color[xpehh.2.ruty.foc$XPEHH_rustica_tytleri > 2.477378] <- '#B03160'
xpehh.2.rugu.foc$color <- 'lightgrey'
xpehh.2.rugu.foc$color[xpehh.2.rugu.foc$XPEHH_rustica_gutturalis < -2.570383] <- '#6BA5CC'
xpehh.2.rugu.foc$color[xpehh.2.rugu.foc$XPEHH_rustica_gutturalis > 2.417505] <- '#B03160'
xpehh.2.guty.foc$color <- 'lightgrey'
xpehh.2.guty.foc$color[xpehh.2.guty.foc$XPEHH_tytleri_gutturalis < -2.542992] <- '#6BA5CC'
xpehh.2.guty.foc$color[xpehh.2.guty.foc$XPEHH_tytleri_gutturalis > 2.447439] <- '#EBB320'

par(mfrow=c(6,1))
## PBS scan
y <- 0.50
y2 <- y-0.1
plot(pbs.50kb.2$fst.50kb.2.ruty.window_pos_1+(pbs.50kb.2$fst.50kb.2.ruty.window_pos_2-pbs.50kb.2$fst.50kb.2.ruty.window_pos_1),pbs.50kb.2$pbs.50kb.2.ru,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="PBS",ylim=c(0,0.6))
lines(pbs.50kb.2$fst.50kb.2.ruty.window_pos_1+(pbs.50kb.2$fst.50kb.2.ruty.window_pos_2-pbs.50kb.2$fst.50kb.2.ruty.window_pos_1),pbs.50kb.2$pbs.50kb.2.ty,lwd=1.5,col='#EBB320',xlim=x)
lines(pbs.50kb.2$fst.50kb.2.ruty.window_pos_1+(pbs.50kb.2$fst.50kb.2.ruty.window_pos_2-pbs.50kb.2$fst.50kb.2.ruty.window_pos_1),pbs.50kb.2$pbs.50kb.2.gu,lwd=1.5,col='#6BA5CC',xlim=x)
for (exon in exons.for.foc){rect(xleft=exons.for.foc$start,ytop=y+0.025,xright=exons.for.foc$end,ybottom=y-0.025,col=alpha('grey',0.75),border=NA)}
for (exon in exons.rev.foc){rect(xleft=exons.rev.foc$start,ytop=y2+0.025,xright=exons.rev.foc$end,ybottom=y2-0.025,col=alpha('grey',0.75),border=NA)}
for (gene in genes.for.foc){arrows(x0=genes.for.foc$start,y0=y,x1=genes.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey50')}
for (gene in genes.rev.foc){arrows(x0=genes.rev.foc$start,y0=y2,x1=genes.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey50')}
for (gene in cand.ts.for.foc){arrows(x0=cand.ts.for.foc$start,y0=y,x1=cand.ts.for.foc$end,y1=y,length=0.035,angle=30,code=2,lwd=1.5,col='grey8')
  text(cand.ts.for.foc$end,y,cand.ts.for.foc$gene,cex=0.65,pos=3)}
for (gene in cand.ts.rev.foc){arrows(x0=cand.ts.rev.foc$start,y0=y2,x1=cand.ts.rev.foc$end,y1=y2,length=0.035,angle=30,code=1,lwd=1.5,col='grey8')
  text(cand.ts.rev.foc$end,y2,cand.ts.rev.foc$gene,cex=0.65,pos=3,)}
## π scan
plot(pi.50kb.2.ru$window_pos_1+(pi.50kb.2.ru$window_pos_2 - pi.50kb.2.ru$window_pos_1)/2,pi.50kb.2.ru$avg_pi,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="π",ylim=c(0,0.004))
lines(pi.50kb.2.ty$window_pos_1+(pi.50kb.2.ty$window_pos_2 - pi.50kb.2.ty$window_pos_1)/2,pi.50kb.2.ty$avg_pi,lwd=1.5,col='#EBB320',xlim=x)
lines(pi.50kb.2.gu$window_pos_1+(pi.50kb.2.gu$window_pos_2 - pi.50kb.2.gu$window_pos_1)/2,pi.50kb.2.gu$avg_pi,lwd=1.5,col='#6BA5CC',xlim=x)
## Tajima's D scan
plot(taj.50kb.2.ru$BIN_START+(taj.50kb.2.ru$BIN_END - taj.50kb.2.ru$BIN_START)/2,taj.50kb.2.ru$TajimaD,type='l',lwd=2,col='#B03160',xlim=x,xlab="",ylab="Tajima's D",ylim=c(-2.5,1.5))
lines(taj.50kb.2.ty$BIN_START+(taj.50kb.2.ty$BIN_END - taj.50kb.2.ty$BIN_START)/2,taj.50kb.2.ty$TajimaD,lwd=1.5,col='#EBB320',xlim=x)
lines(taj.50kb.2.gu$BIN_START+(taj.50kb.2.gu$BIN_END - taj.50kb.2.gu$BIN_START)/2,taj.50kb.2.gu$TajimaD,lwd=1.5,col='#6BA5CC',xlim=x)
abline(h=mean(taj.10kb.ru$TajimaD),lty=2,col=alpha('#B03160',0.5))
abline(h=mean(taj.10kb.ty$TajimaD),lty=2,col=alpha('#EBB320',0.5))
abline(h=mean(taj.10kb.gu$TajimaD),lty=2,col=alpha('#6BA5CC',0.5))
## xp-EHH scans
plot(xpehh.2.ruty.foc$POSITION,xpehh.2.ruty.foc$XPEHH_rustica_tytleri,pch=20,cex=0.5,col=xpehh.2.ruty.foc$color,xlim=x,ylim=c(-7,7),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.2.rugu.foc$POSITION,xpehh.2.rugu.foc$XPEHH_rustica_gutturalis,pch=20,cex=0.5,col=xpehh.2.rugu.foc$color,xlim=x,ylim=c(-7,7),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')
plot(xpehh.2.guty.foc$POSITION,xpehh.2.guty.foc$XPEHH_tytleri_gutturalis,pch=20,cex=0.5,col=xpehh.2.guty.foc$color,xlim=x,ylim=c(-7,7),ylab='xp-EHH')
abline(h=0,lty=2,col='grey2')

