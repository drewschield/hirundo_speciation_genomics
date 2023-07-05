############################################################################
# Barn swallow genome-wide association mapping - LMM results
############################################################################

# We ran genotype x phenotype association mapping analysis in GEMMA to find
# associations between genomic regions and trait variation.

# This script contains commands for parsing, plotting, and performing summary
# statistics on LMM results.

# We'll examine the results of analyses including individuals from hybrid
# zones.

# Conventions: 'vc' = ventral color; 'ts' = tail streamer

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/gemma/')

library(data.table)
library(tidyverse)
library(dplyr)
library(scales)
library(Rmisc)
library(matrixStats)
library(qqman)

### Read in and format genome-wide LMM results-------------------------------

# Hybrid results
vc.hy <- read.table('./lmm_full_prune/gwas_full_lmm.hybrid-breast-bright.assoc.prune.txt',header=T)
ts.hy <- read.table('./lmm_full_prune/gwas_full_lmm.hybrid-tail-streamer.assoc.prune.txt',header=T)

# 'Full' results
vc.full <- read.table('./lmm_full_prune/gwas_full_lmm.full-breast-bright.assoc.prune.txt',header=T)
ts.full <- read.table('./lmm_full_prune/gwas_full_lmm.full-tail-streamer.assoc.prune.txt',header=T)

# Males only results
vc.male <- read.table('./lmm_full_prune/gwas_full_lmm.male-breast-bright.assoc.prune.txt',header=T)
ts.male <- read.table('./lmm_full_prune/gwas_full_lmm.male-tail-streamer.assoc.prune.txt',header=T)

# Females only results
vc.female <- read.table('./lmm_full_prune/gwas_full_lmm.female-breast-bright.assoc.prune.txt',header=T)
ts.female <- read.table('./lmm_full_prune/gwas_full_lmm.female-tail-streamer.assoc.prune.txt',header=T)

# 'Full' results with randomized phenotypes
vc.rand <- read.table('./lmm_full_prune/gwas_full_lmm.full-breast-bright-random.assoc.prune.txt',header=T)
ts.rand <- read.table('./lmm_full_prune/gwas_full_lmm.full-tail-streamer-random.assoc.prune.txt',header=T)

# load file to rename chromosome names
chr_rename <- read.table('./chr_rename.txt',header=T)

## Reformat with integer chromosome labels

vc.hy <- bind_cols(vc.hy) %>%
  rename(CHR = chr)
vc.hy2 <- left_join(vc.hy, chr_rename, by="CHR")

ts.hy <- bind_cols(ts.hy) %>%
  rename(CHR = chr)
ts.hy2 <- left_join(ts.hy, chr_rename, by="CHR")

vc.full <- bind_cols(vc.full) %>%
  rename(CHR = chr)
vc.full2 <- left_join(vc.full, chr_rename, by="CHR")

ts.full <- bind_cols(ts.full) %>%
  rename(CHR = chr)
ts.full2 <- left_join(ts.full, chr_rename, by="CHR")

vc.male <- bind_cols(vc.male) %>%
  rename(CHR = chr)
vc.male2 <- left_join(vc.male, chr_rename, by="CHR")

ts.male <- bind_cols(ts.male) %>%
  rename(CHR = chr)
ts.male2 <- left_join(ts.male, chr_rename, by="CHR")

vc.female <- bind_cols(vc.female) %>%
  rename(CHR = chr)
vc.female2 <- left_join(vc.female, chr_rename, by="CHR")

ts.female <- bind_cols(ts.female) %>%
  rename(CHR = chr)
ts.female2 <- left_join(ts.female, chr_rename, by="CHR")

vc.rand <- bind_cols(vc.rand) %>%
  rename(CHR = chr)
vc.rand2 <- left_join(vc.rand, chr_rename, by="CHR")

ts.rand <- bind_cols(ts.rand) %>%
  rename(CHR = chr)
ts.rand2 <- left_join(ts.rand, chr_rename, by="CHR")

## Calculate and append -log10 Wald p-value

# Ventral color
logP <- -log10(vc.hy2$p_wald)
vc.hy2$logP <- logP

logP <- -log10(ts.hy2$p_wald)
ts.hy2$logP <- logP

logP <- -log10(vc.full2$p_wald)
vc.full2$logP <- logP

logP <- -log10(ts.full2$p_wald)
ts.full2$logP <- logP

logP <- -log10(vc.male2$p_wald)
vc.male2$logP <- logP

logP <- -log10(ts.male2$p_wald)
ts.male2$logP <- logP

logP <- -log10(vc.female2$p_wald)
vc.female2$logP <- logP

logP <- -log10(ts.female2$p_wald)
ts.female2$logP <- logP

logP <- -log10(vc.rand2$p_wald)
vc.rand2$logP <- logP

logP <- -log10(ts.rand2$p_wald)
ts.rand2$logP <- logP

### Genome-wide Manhattan plots---------------------------------------------

# Numbers of SNPs in analyses (for Bonferroni correction):

# Ventral color
# hybrid  = 9033285
# full    = 9311933 
# male    = 8851674
# female  = 8758246
# random  = 9311933

# Tail streamer
# hybrid  = 8870504
# full    = 9246918 
# male    = 8744006
# female  = 8773709
# random  = 9246918

colors <- c("#000000", "#7a7a7a", "#adadad")

# Output individually at 4 x 10
par(mfrow=c(1,1))
manhattan(filter(vc.hy2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/9033285), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(vc.full2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/9311933), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(vc.male2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/8851674), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(vc.female2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/8758246), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(vc.rand2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/9311933), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))

manhattan(filter(ts.hy2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/8870504), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(ts.full2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/9246918), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(ts.male2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/8744006), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(ts.female2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/8773709), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))
manhattan(filter(ts.rand2,logP>2), chr="CHR_NUM", snp = "ps", bp="ps", p="logP", logp=FALSE, genomewideline=FALSE, suggestiveline=-log10(0.05/9246918), ylab="-log(P)", xlab="Chromosome", col=alpha(colors,0.25),ylim=c(0,20))

### Plot chromosome-specific results-----------------------------------------

# Read in data for focal chromosomes
vc.hy.chr1a <- read.table('./lmm_full_chrom/gwas_full_lmm.hybrid-breast-bright.assoc.chr1A-NC_053453.1.txt',header=T)
vc.hy.chrz <- read.table('./lmm_full_chrom/gwas_full_lmm.hybrid-breast-bright.assoc.chrZ-NC_053488.1.txt',header=T)
ts.hy.chr2 <- read.table('./lmm_full_chrom/gwas_full_lmm.hybrid-tail-streamer.assoc.chr2-NC_053450.1.txt',header=T)
# Get rid of super low stuff for plotting's sake
vc.hy.chr1a <- vc.hy.chr1a[which(-log10(vc.hy.chr1a$p_wald)>=1),]
vc.hy.chrz <- vc.hy.chrz[which(-log10(vc.hy.chrz$p_wald)>=1),]
ts.hy.chr2 <- ts.hy.chr2[which(-log10(ts.hy.chr2$p_wald)>=1),]

## Plot LMM results

# Ventral color
# Set colors below/above bonferroni-corrected threshold
vc.hy.chr1a$color <- 'black'
vc.hy.chr1a$color[vc.hy.chr1a$p_wald<=0.05/9033285] <- 'brown'
vc.hy.chr1a$color[vc.hy.chr1a$p_wald>0.05/9033285] <- 'lightgrey'

vc.hy.chrz$color <- 'black'
vc.hy.chrz$color[vc.hy.chrz$p_wald<=0.05/9033285] <- 'brown'
vc.hy.chrz$color[vc.hy.chrz$p_wald>0.05/9033285] <- 'lightgrey'

# Plot LMM panels
# For ventral color, we'll set Chr1A and Z to common chromosome 2 x-axis for scale
par(mfrow=c(1,1))
# Output at 4 x 8
plot(vc.hy.chr1a$ps,-log10(vc.hy.chr1a$p_wald),pch=20,cex=0.6,col=vc.hy.chr1a$color,ylab='-log10(P)',xlab='Chromosome Position (Mb)',xlim=c(0,156035725))
plot(vc.hy.chrz$ps,-log10(vc.hy.chrz$p_wald),pch=20,cex=0.6,col=vc.hy.chrz$color,ylab='-log10(P)',xlab='Chromosome Position (Mb)',xlim=c(0,156035725),ylim=c(1,20))

# Tail streamer
# Set colors below/above bonferroni-corrected threshold
ts.hy.chr2$color <- 'black'
ts.hy.chr2$color[ts.hy.chr2$p_wald<=0.05/8870504] <- 'darkblue'
ts.hy.chr2$color[ts.hy.chr2$p_wald>0.05/8870504] <- 'lightgrey'

# Plot LMM panels
# For ventral color, we'll set Chr1A and Z to common x-axis for scale
par(mfrow=c(1,1))
# Output at 4 x 8
plot(ts.hy.chr2$ps,-log10(ts.hy.chr2$p_wald),pch=20,cex=0.6,col=ts.hy.chr2$color,ylab='-log10(P)',xlab='Chromosome Position (Mb)')
