############################################################################
# Barn swallow genome-wide association mapping - BSLMM results
############################################################################

# We ran genotype x phenotype association mapping analysis in GEMMA to find
# associations between genomic regions and trait variation.

# This script contains commands for parsing, plotting, and performing summary
# statistics on BSLMM results.

# The first section involves plotting and summarization of hyperparameters
# that estimate the quantitative genetics of the traits, including the total
# proportion of phenotypic variation explained by all SNPs in our model, and
# the posterior inclusion probability (PIP) that specific SNPs have non-zero
# effects on phenotype.

# We'll examine the results of analyses including individuals from hybrid
# zones.

# Conventions: 'vc' = ventral color; 'ts' = tail streamer

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/gemma/bslmm/')

library(data.table)
library(tidyverse)
library(dplyr)
library(scales)
library(Rmisc)
library(matrixStats)
library(qqman)

### Read in hyperparameter data---------------------------------------------

hyp.vc <- read.table('./gwas_full_bslmm.hybrid-breast-bright.run_hyp.txt',header=T)
hyp.ts <- read.table('./gwas_full_bslmm.hybrid-tail-streamer.run_hyp.txt',header=T)

### Read in parameter data--------------------------------------------------

par.vc <- read.table('./gwas_full_bslmm.hybrid-breast-bright.run_sparse.txt',header=T)
par.vc <- par.vc[which(par.vc$gamma>=0.01),]
par.ts <- read.table('./gwas_full_bslmm.hybrid-tail-streamer.run_sparse.txt',header=T)
par.ts <- par.ts[which(par.ts$gamma>=0.01),]

### Plot hyperparameter posterior distributions-----------------------------

## Output at 3 x 8

par(mfrow=c(1,3))
plot(density(hyp.vc$pve),lwd=2,xlab='PVE',main=NA)
plot(density(hyp.vc$pge),lwd=2,xlab='PGE',main=NA)
plot(density(hyp.vc$n_gamma),lwd=2,xlab='n SNPs',main=NA)

plot(density(hyp.ts$pve),lwd=2,xlab='PVE',main=NA)
plot(density(hyp.ts$pge),lwd=2,xlab='PGE',main=NA)
plot(density(hyp.ts$n_gamma,bw=1),lwd=2,xlab='n SNPs',main=NA)

### Parse parameters for SNPs with sparse effects (ß > 0)-------------------

# Note: this is gross and should be put in loops.

## Ventral color
vc.chr1 <- par.vc[which(par.vc$CHR=='NC_053451.1'),]
vc.chr1a <- par.vc[which(par.vc$CHR=='NC_053453.1'),]
vc.chr2 <- par.vc[which(par.vc$CHR=='NC_053450.1'),]
vc.chr3 <- par.vc[which(par.vc$CHR=='NC_053452.1'),]
vc.chr4 <- par.vc[which(par.vc$CHR=='NC_053454.1'),]
vc.chr4a <- par.vc[which(par.vc$CHR=='NC_053470.1'),]
vc.chr5 <- par.vc[which(par.vc$CHR=='NC_053455.1'),]
vc.chr6 <- par.vc[which(par.vc$CHR=='NC_053457.1'),]
vc.chr7 <- par.vc[which(par.vc$CHR=='NC_053456.1'),]
vc.chr8 <- par.vc[which(par.vc$CHR=='NC_053458.1'),]
vc.chr9 <- par.vc[which(par.vc$CHR=='NC_053459.1'),]
vc.chr10 <- par.vc[which(par.vc$CHR=='NC_053462.1'),]
vc.chr11 <- par.vc[which(par.vc$CHR=='NC_053460.1'),]
vc.chr12 <- par.vc[which(par.vc$CHR=='NC_053461.1'),]
vc.chr13 <- par.vc[which(par.vc$CHR=='NC_053463.1'),]
vc.chr14 <- par.vc[which(par.vc$CHR=='NC_053464.1'),]
vc.chr15 <- par.vc[which(par.vc$CHR=='NC_053466.1'),]
vc.chr17 <- par.vc[which(par.vc$CHR=='NC_053469.1'),]
vc.chr18 <- par.vc[which(par.vc$CHR=='NC_053467.1'),]
vc.chr19 <- par.vc[which(par.vc$CHR=='NC_053468.1'),]
vc.chr20 <- par.vc[which(par.vc$CHR=='NC_053465.1'),]
vc.chr21 <- par.vc[which(par.vc$CHR=='NC_053471.1'),]
vc.chr22 <- par.vc[which(par.vc$CHR=='NC_053477.1'),]
vc.chr23 <- par.vc[which(par.vc$CHR=='NC_053474.1'),]
vc.chr24 <- par.vc[which(par.vc$CHR=='NC_053472.1'),]
vc.chr25 <- par.vc[which(par.vc$CHR=='NC_053478.1'),]
vc.chr26 <- par.vc[which(par.vc$CHR=='NC_053473.1'),]
vc.chr27 <- par.vc[which(par.vc$CHR=='NC_053476.1'),]
vc.chr28 <- par.vc[which(par.vc$CHR=='NC_053475.1'),]
vc.chr29 <- par.vc[which(par.vc$CHR=='NC_053479.1'),]
vc.chr31 <- par.vc[which(par.vc$CHR=='NC_053480.1'),]
vc.chr32 <- par.vc[which(par.vc$CHR=='NC_053481.1'),]
vc.chr33 <- par.vc[which(par.vc$CHR=='NC_053482.1'),]
vc.chr34 <- par.vc[which(par.vc$CHR=='NC_053483.1'),]
vc.chr35 <- par.vc[which(par.vc$CHR=='NC_053484.1'),]
vc.chr36 <- par.vc[which(par.vc$CHR=='NC_053485.1'),]
vc.chr37 <- par.vc[which(par.vc$CHR=='NC_053486.1'),]
vc.chrz <- par.vc[which(par.vc$CHR=='NC_053488.1'),]

vc.means <- c()
vc.means=append(vc.means, mean(vc.chr1$gamma))
vc.means=append(vc.means, mean(vc.chr1a$gamma))
vc.means=append(vc.means, mean(vc.chr2$gamma))
vc.means=append(vc.means, mean(vc.chr3$gamma))
vc.means=append(vc.means, mean(vc.chr4$gamma))
vc.means=append(vc.means, mean(vc.chr4a$gamma))
vc.means=append(vc.means, mean(vc.chr5$gamma))
vc.means=append(vc.means, mean(vc.chr6$gamma))
vc.means=append(vc.means, mean(vc.chr7$gamma))
vc.means=append(vc.means, mean(vc.chr8$gamma))
vc.means=append(vc.means, mean(vc.chr9$gamma))
vc.means=append(vc.means, mean(vc.chr10$gamma))
vc.means=append(vc.means, mean(vc.chr11$gamma))
vc.means=append(vc.means, mean(vc.chr12$gamma))
vc.means=append(vc.means, mean(vc.chr13$gamma))
vc.means=append(vc.means, mean(vc.chr14$gamma))
vc.means=append(vc.means, mean(vc.chr15$gamma))
vc.means=append(vc.means, mean(vc.chr17$gamma))
vc.means=append(vc.means, mean(vc.chr18$gamma))
vc.means=append(vc.means, mean(vc.chr19$gamma))
vc.means=append(vc.means, mean(vc.chr20$gamma))
vc.means=append(vc.means, mean(vc.chr21$gamma))
vc.means=append(vc.means, mean(vc.chr22$gamma))
vc.means=append(vc.means, mean(vc.chr23$gamma))
vc.means=append(vc.means, mean(vc.chr24$gamma))
vc.means=append(vc.means, mean(vc.chr25$gamma))
vc.means=append(vc.means, mean(vc.chr26$gamma))
vc.means=append(vc.means, mean(vc.chr27$gamma))
vc.means=append(vc.means, mean(vc.chr28$gamma))
vc.means=append(vc.means, mean(vc.chr29$gamma))
vc.means=append(vc.means, mean(vc.chr31$gamma))
vc.means=append(vc.means, mean(vc.chr32$gamma))
vc.means=append(vc.means, mean(vc.chr33$gamma))
vc.means=append(vc.means, mean(vc.chr34$gamma))
vc.means=append(vc.means, mean(vc.chr35$gamma))
vc.means=append(vc.means, mean(vc.chr36$gamma))
vc.means=append(vc.means, mean(vc.chr37$gamma))
vc.means=append(vc.means, mean(vc.chrz$gamma))

std <- function(x) sd(x)/sqrt(length(x))

vc.stds <- c()
vc.stds=append(vc.stds,std(vc.chr1$gamma))
vc.stds=append(vc.stds,std(vc.chr1a$gamma))
vc.stds=append(vc.stds,std(vc.chr2$gamma))
vc.stds=append(vc.stds,std(vc.chr3$gamma))
vc.stds=append(vc.stds,std(vc.chr4$gamma))
vc.stds=append(vc.stds,std(vc.chr4a$gamma))
vc.stds=append(vc.stds,std(vc.chr5$gamma))
vc.stds=append(vc.stds,std(vc.chr6$gamma))
vc.stds=append(vc.stds,std(vc.chr7$gamma))
vc.stds=append(vc.stds,std(vc.chr8$gamma))
vc.stds=append(vc.stds,std(vc.chr9$gamma))
vc.stds=append(vc.stds,std(vc.chr10$gamma))
vc.stds=append(vc.stds,std(vc.chr11$gamma))
vc.stds=append(vc.stds,std(vc.chr12$gamma))
vc.stds=append(vc.stds,std(vc.chr13$gamma))
vc.stds=append(vc.stds,std(vc.chr14$gamma))
vc.stds=append(vc.stds,std(vc.chr15$gamma))
vc.stds=append(vc.stds,std(vc.chr17$gamma))
vc.stds=append(vc.stds,std(vc.chr18$gamma))
vc.stds=append(vc.stds,std(vc.chr19$gamma))
vc.stds=append(vc.stds,std(vc.chr20$gamma))
vc.stds=append(vc.stds,std(vc.chr21$gamma))
vc.stds=append(vc.stds,std(vc.chr22$gamma))
vc.stds=append(vc.stds,std(vc.chr23$gamma))
vc.stds=append(vc.stds,std(vc.chr24$gamma))
vc.stds=append(vc.stds,std(vc.chr25$gamma))
vc.stds=append(vc.stds,std(vc.chr26$gamma))
vc.stds=append(vc.stds,std(vc.chr27$gamma))
vc.stds=append(vc.stds,std(vc.chr28$gamma))
vc.stds=append(vc.stds,std(vc.chr29$gamma))
vc.stds=append(vc.stds,std(vc.chr31$gamma))
vc.stds=append(vc.stds,std(vc.chr32$gamma))
vc.stds=append(vc.stds,std(vc.chr33$gamma))
vc.stds=append(vc.stds,std(vc.chr34$gamma))
vc.stds=append(vc.stds,std(vc.chr35$gamma))
vc.stds=append(vc.stds,std(vc.chr36$gamma))
vc.stds=append(vc.stds,std(vc.chr37$gamma))
vc.stds=append(vc.stds,std(vc.chrz$gamma))

vc.maxs <- c()
vc.maxs=append(vc.maxs,max(vc.chr1$gamma))
vc.maxs=append(vc.maxs,max(vc.chr1a$gamma))
vc.maxs=append(vc.maxs,max(vc.chr2$gamma))
vc.maxs=append(vc.maxs,max(vc.chr3$gamma))
vc.maxs=append(vc.maxs,max(vc.chr4$gamma))
vc.maxs=append(vc.maxs,max(vc.chr4a$gamma))
vc.maxs=append(vc.maxs,max(vc.chr5$gamma))
vc.maxs=append(vc.maxs,max(vc.chr6$gamma))
vc.maxs=append(vc.maxs,max(vc.chr7$gamma))
vc.maxs=append(vc.maxs,max(vc.chr8$gamma))
vc.maxs=append(vc.maxs,max(vc.chr9$gamma))
vc.maxs=append(vc.maxs,max(vc.chr10$gamma))
vc.maxs=append(vc.maxs,max(vc.chr11$gamma))
vc.maxs=append(vc.maxs,max(vc.chr12$gamma))
vc.maxs=append(vc.maxs,max(vc.chr13$gamma))
vc.maxs=append(vc.maxs,max(vc.chr14$gamma))
vc.maxs=append(vc.maxs,max(vc.chr15$gamma))
vc.maxs=append(vc.maxs,max(vc.chr17$gamma))
vc.maxs=append(vc.maxs,max(vc.chr18$gamma))
vc.maxs=append(vc.maxs,max(vc.chr19$gamma))
vc.maxs=append(vc.maxs,max(vc.chr20$gamma))
vc.maxs=append(vc.maxs,max(vc.chr21$gamma))
vc.maxs=append(vc.maxs,max(vc.chr22$gamma))
vc.maxs=append(vc.maxs,max(vc.chr23$gamma))
vc.maxs=append(vc.maxs,max(vc.chr24$gamma))
vc.maxs=append(vc.maxs,max(vc.chr25$gamma))
vc.maxs=append(vc.maxs,max(vc.chr26$gamma))
vc.maxs=append(vc.maxs,max(vc.chr27$gamma))
vc.maxs=append(vc.maxs,max(vc.chr28$gamma))
vc.maxs=append(vc.maxs,max(vc.chr29$gamma))
vc.maxs=append(vc.maxs,max(vc.chr31$gamma))
vc.maxs=append(vc.maxs,max(vc.chr32$gamma))
vc.maxs=append(vc.maxs,max(vc.chr33$gamma))
vc.maxs=append(vc.maxs,max(vc.chr34$gamma))
vc.maxs=append(vc.maxs,max(vc.chr35$gamma))
vc.maxs=append(vc.maxs,max(vc.chr36$gamma))
vc.maxs=append(vc.maxs,max(vc.chr37$gamma))
vc.maxs=append(vc.maxs,max(vc.chrz$gamma))

vc.sums <- c()
vc.sums=append(vc.sums,sum(vc.chr1$gamma))
vc.sums=append(vc.sums,sum(vc.chr1a$gamma))
vc.sums=append(vc.sums,sum(vc.chr2$gamma))
vc.sums=append(vc.sums,sum(vc.chr3$gamma))
vc.sums=append(vc.sums,sum(vc.chr4$gamma))
vc.sums=append(vc.sums,sum(vc.chr4a$gamma))
vc.sums=append(vc.sums,sum(vc.chr5$gamma))
vc.sums=append(vc.sums,sum(vc.chr6$gamma))
vc.sums=append(vc.sums,sum(vc.chr7$gamma))
vc.sums=append(vc.sums,sum(vc.chr8$gamma))
vc.sums=append(vc.sums,sum(vc.chr9$gamma))
vc.sums=append(vc.sums,sum(vc.chr10$gamma))
vc.sums=append(vc.sums,sum(vc.chr11$gamma))
vc.sums=append(vc.sums,sum(vc.chr12$gamma))
vc.sums=append(vc.sums,sum(vc.chr13$gamma))
vc.sums=append(vc.sums,sum(vc.chr14$gamma))
vc.sums=append(vc.sums,sum(vc.chr15$gamma))
vc.sums=append(vc.sums,sum(vc.chr17$gamma))
vc.sums=append(vc.sums,sum(vc.chr18$gamma))
vc.sums=append(vc.sums,sum(vc.chr19$gamma))
vc.sums=append(vc.sums,sum(vc.chr20$gamma))
vc.sums=append(vc.sums,sum(vc.chr21$gamma))
vc.sums=append(vc.sums,sum(vc.chr22$gamma))
vc.sums=append(vc.sums,sum(vc.chr23$gamma))
vc.sums=append(vc.sums,sum(vc.chr24$gamma))
vc.sums=append(vc.sums,sum(vc.chr25$gamma))
vc.sums=append(vc.sums,sum(vc.chr26$gamma))
vc.sums=append(vc.sums,sum(vc.chr27$gamma))
vc.sums=append(vc.sums,sum(vc.chr28$gamma))
vc.sums=append(vc.sums,sum(vc.chr29$gamma))
vc.sums=append(vc.sums,sum(vc.chr31$gamma))
vc.sums=append(vc.sums,sum(vc.chr32$gamma))
vc.sums=append(vc.sums,sum(vc.chr33$gamma))
vc.sums=append(vc.sums,sum(vc.chr34$gamma))
vc.sums=append(vc.sums,sum(vc.chr35$gamma))
vc.sums=append(vc.sums,sum(vc.chr36$gamma))
vc.sums=append(vc.sums,sum(vc.chr37$gamma))
vc.sums=append(vc.sums,sum(vc.chrz$gamma))

## Tail streamer
ts.chr1 <- par.ts[which(par.ts$CHR=='NC_053451.1'),]
ts.chr1a <- par.ts[which(par.ts$CHR=='NC_053453.1'),]
ts.chr2 <- par.ts[which(par.ts$CHR=='NC_053450.1'),]
ts.chr3 <- par.ts[which(par.ts$CHR=='NC_053452.1'),]
ts.chr4 <- par.ts[which(par.ts$CHR=='NC_053454.1'),]
ts.chr4a <- par.ts[which(par.ts$CHR=='NC_053470.1'),]
ts.chr5 <- par.ts[which(par.ts$CHR=='NC_053455.1'),]
ts.chr6 <- par.ts[which(par.ts$CHR=='NC_053457.1'),]
ts.chr7 <- par.ts[which(par.ts$CHR=='NC_053456.1'),]
ts.chr8 <- par.ts[which(par.ts$CHR=='NC_053458.1'),]
ts.chr9 <- par.ts[which(par.ts$CHR=='NC_053459.1'),]
ts.chr10 <- par.ts[which(par.ts$CHR=='NC_053462.1'),]
ts.chr11 <- par.ts[which(par.ts$CHR=='NC_053460.1'),]
ts.chr12 <- par.ts[which(par.ts$CHR=='NC_053461.1'),]
ts.chr13 <- par.ts[which(par.ts$CHR=='NC_053463.1'),]
ts.chr14 <- par.ts[which(par.ts$CHR=='NC_053464.1'),]
ts.chr15 <- par.ts[which(par.ts$CHR=='NC_053466.1'),]
ts.chr17 <- par.ts[which(par.ts$CHR=='NC_053469.1'),]
ts.chr18 <- par.ts[which(par.ts$CHR=='NC_053467.1'),]
ts.chr19 <- par.ts[which(par.ts$CHR=='NC_053468.1'),]
ts.chr20 <- par.ts[which(par.ts$CHR=='NC_053465.1'),]
ts.chr21 <- par.ts[which(par.ts$CHR=='NC_053471.1'),]
ts.chr22 <- par.ts[which(par.ts$CHR=='NC_053477.1'),]
ts.chr23 <- par.ts[which(par.ts$CHR=='NC_053474.1'),]
ts.chr24 <- par.ts[which(par.ts$CHR=='NC_053472.1'),]
ts.chr25 <- par.ts[which(par.ts$CHR=='NC_053478.1'),]
ts.chr26 <- par.ts[which(par.ts$CHR=='NC_053473.1'),]
ts.chr27 <- par.ts[which(par.ts$CHR=='NC_053476.1'),]
ts.chr28 <- par.ts[which(par.ts$CHR=='NC_053475.1'),]
ts.chr29 <- par.ts[which(par.ts$CHR=='NC_053479.1'),]
ts.chr31 <- par.ts[which(par.ts$CHR=='NC_053480.1'),]
ts.chr32 <- par.ts[which(par.ts$CHR=='NC_053481.1'),]
ts.chr33 <- par.ts[which(par.ts$CHR=='NC_053482.1'),]
ts.chr34 <- par.ts[which(par.ts$CHR=='NC_053483.1'),]
ts.chr35 <- par.ts[which(par.ts$CHR=='NC_053484.1'),]
ts.chr36 <- par.ts[which(par.ts$CHR=='NC_053485.1'),]
ts.chr37 <- par.ts[which(par.ts$CHR=='NC_053486.1'),]
ts.chrz <- par.ts[which(par.ts$CHR=='NC_053488.1'),]

ts.means <- c()
ts.means=append(ts.means, mean(ts.chr1$gamma))
ts.means=append(ts.means, mean(ts.chr1a$gamma))
ts.means=append(ts.means, mean(ts.chr2$gamma))
ts.means=append(ts.means, mean(ts.chr3$gamma))
ts.means=append(ts.means, mean(ts.chr4$gamma))
ts.means=append(ts.means, mean(ts.chr4a$gamma))
ts.means=append(ts.means, mean(ts.chr5$gamma))
ts.means=append(ts.means, mean(ts.chr6$gamma))
ts.means=append(ts.means, mean(ts.chr7$gamma))
ts.means=append(ts.means, mean(ts.chr8$gamma))
ts.means=append(ts.means, mean(ts.chr9$gamma))
ts.means=append(ts.means, mean(ts.chr10$gamma))
ts.means=append(ts.means, mean(ts.chr11$gamma))
ts.means=append(ts.means, mean(ts.chr12$gamma))
ts.means=append(ts.means, mean(ts.chr13$gamma))
ts.means=append(ts.means, mean(ts.chr14$gamma))
ts.means=append(ts.means, mean(ts.chr15$gamma))
ts.means=append(ts.means, mean(ts.chr17$gamma))
ts.means=append(ts.means, mean(ts.chr18$gamma))
ts.means=append(ts.means, mean(ts.chr19$gamma))
ts.means=append(ts.means, mean(ts.chr20$gamma))
ts.means=append(ts.means, mean(ts.chr21$gamma))
ts.means=append(ts.means, mean(ts.chr22$gamma))
ts.means=append(ts.means, mean(ts.chr23$gamma))
ts.means=append(ts.means, mean(ts.chr24$gamma))
ts.means=append(ts.means, mean(ts.chr25$gamma))
ts.means=append(ts.means, mean(ts.chr26$gamma))
ts.means=append(ts.means, mean(ts.chr27$gamma))
ts.means=append(ts.means, mean(ts.chr28$gamma))
ts.means=append(ts.means, mean(ts.chr29$gamma))
ts.means=append(ts.means, mean(ts.chr31$gamma))
ts.means=append(ts.means, mean(ts.chr32$gamma))
ts.means=append(ts.means, mean(ts.chr33$gamma))
ts.means=append(ts.means, mean(ts.chr34$gamma))
ts.means=append(ts.means, mean(ts.chr35$gamma))
ts.means=append(ts.means, mean(ts.chr36$gamma))
ts.means=append(ts.means, mean(ts.chr37$gamma))
ts.means=append(ts.means, mean(ts.chrz$gamma))

std <- function(x) sd(x)/sqrt(length(x))

ts.stds <- c()
ts.stds=append(ts.stds,std(ts.chr1$gamma))
ts.stds=append(ts.stds,std(ts.chr1a$gamma))
ts.stds=append(ts.stds,std(ts.chr2$gamma))
ts.stds=append(ts.stds,std(ts.chr3$gamma))
ts.stds=append(ts.stds,std(ts.chr4$gamma))
ts.stds=append(ts.stds,std(ts.chr4a$gamma))
ts.stds=append(ts.stds,std(ts.chr5$gamma))
ts.stds=append(ts.stds,std(ts.chr6$gamma))
ts.stds=append(ts.stds,std(ts.chr7$gamma))
ts.stds=append(ts.stds,std(ts.chr8$gamma))
ts.stds=append(ts.stds,std(ts.chr9$gamma))
ts.stds=append(ts.stds,std(ts.chr10$gamma))
ts.stds=append(ts.stds,std(ts.chr11$gamma))
ts.stds=append(ts.stds,std(ts.chr12$gamma))
ts.stds=append(ts.stds,std(ts.chr13$gamma))
ts.stds=append(ts.stds,std(ts.chr14$gamma))
ts.stds=append(ts.stds,std(ts.chr15$gamma))
ts.stds=append(ts.stds,std(ts.chr17$gamma))
ts.stds=append(ts.stds,std(ts.chr18$gamma))
ts.stds=append(ts.stds,std(ts.chr19$gamma))
ts.stds=append(ts.stds,std(ts.chr20$gamma))
ts.stds=append(ts.stds,std(ts.chr21$gamma))
ts.stds=append(ts.stds,std(ts.chr22$gamma))
ts.stds=append(ts.stds,std(ts.chr23$gamma))
ts.stds=append(ts.stds,std(ts.chr24$gamma))
ts.stds=append(ts.stds,std(ts.chr25$gamma))
ts.stds=append(ts.stds,std(ts.chr26$gamma))
ts.stds=append(ts.stds,std(ts.chr27$gamma))
ts.stds=append(ts.stds,std(ts.chr28$gamma))
ts.stds=append(ts.stds,std(ts.chr29$gamma))
ts.stds=append(ts.stds,std(ts.chr31$gamma))
ts.stds=append(ts.stds,std(ts.chr32$gamma))
ts.stds=append(ts.stds,std(ts.chr33$gamma))
ts.stds=append(ts.stds,std(ts.chr34$gamma))
ts.stds=append(ts.stds,std(ts.chr35$gamma))
ts.stds=append(ts.stds,std(ts.chr36$gamma))
ts.stds=append(ts.stds,std(ts.chr37$gamma))
ts.stds=append(ts.stds,std(ts.chrz$gamma))

ts.maxs <- c()
ts.maxs=append(ts.maxs,max(ts.chr1$gamma))
ts.maxs=append(ts.maxs,max(ts.chr1a$gamma))
ts.maxs=append(ts.maxs,max(ts.chr2$gamma))
ts.maxs=append(ts.maxs,max(ts.chr3$gamma))
ts.maxs=append(ts.maxs,max(ts.chr4$gamma))
ts.maxs=append(ts.maxs,max(ts.chr4a$gamma))
ts.maxs=append(ts.maxs,max(ts.chr5$gamma))
ts.maxs=append(ts.maxs,max(ts.chr6$gamma))
ts.maxs=append(ts.maxs,max(ts.chr7$gamma))
ts.maxs=append(ts.maxs,max(ts.chr8$gamma))
ts.maxs=append(ts.maxs,max(ts.chr9$gamma))
ts.maxs=append(ts.maxs,max(ts.chr10$gamma))
ts.maxs=append(ts.maxs,max(ts.chr11$gamma))
ts.maxs=append(ts.maxs,max(ts.chr12$gamma))
ts.maxs=append(ts.maxs,max(ts.chr13$gamma))
ts.maxs=append(ts.maxs,max(ts.chr14$gamma))
ts.maxs=append(ts.maxs,max(ts.chr15$gamma))
ts.maxs=append(ts.maxs,max(ts.chr17$gamma))
ts.maxs=append(ts.maxs,max(ts.chr18$gamma))
ts.maxs=append(ts.maxs,max(ts.chr19$gamma))
ts.maxs=append(ts.maxs,max(ts.chr20$gamma))
ts.maxs=append(ts.maxs,max(ts.chr21$gamma))
ts.maxs=append(ts.maxs,max(ts.chr22$gamma))
ts.maxs=append(ts.maxs,max(ts.chr23$gamma))
ts.maxs=append(ts.maxs,max(ts.chr24$gamma))
ts.maxs=append(ts.maxs,max(ts.chr25$gamma))
ts.maxs=append(ts.maxs,max(ts.chr26$gamma))
ts.maxs=append(ts.maxs,max(ts.chr27$gamma))
ts.maxs=append(ts.maxs,max(ts.chr28$gamma))
ts.maxs=append(ts.maxs,max(ts.chr29$gamma))
ts.maxs=append(ts.maxs,max(ts.chr31$gamma))
ts.maxs=append(ts.maxs,max(ts.chr32$gamma))
ts.maxs=append(ts.maxs,max(ts.chr33$gamma))
ts.maxs=append(ts.maxs,max(ts.chr34$gamma))
ts.maxs=append(ts.maxs,max(ts.chr35$gamma))
ts.maxs=append(ts.maxs,max(ts.chr36$gamma))
ts.maxs=append(ts.maxs,max(ts.chr37$gamma))
ts.maxs=append(ts.maxs,max(ts.chrz$gamma))

ts.sums <- c()
ts.sums=append(ts.sums,sum(ts.chr1$gamma))
ts.sums=append(ts.sums,sum(ts.chr1a$gamma))
ts.sums=append(ts.sums,sum(ts.chr2$gamma))
ts.sums=append(ts.sums,sum(ts.chr3$gamma))
ts.sums=append(ts.sums,sum(ts.chr4$gamma))
ts.sums=append(ts.sums,sum(ts.chr4a$gamma))
ts.sums=append(ts.sums,sum(ts.chr5$gamma))
ts.sums=append(ts.sums,sum(ts.chr6$gamma))
ts.sums=append(ts.sums,sum(ts.chr7$gamma))
ts.sums=append(ts.sums,sum(ts.chr8$gamma))
ts.sums=append(ts.sums,sum(ts.chr9$gamma))
ts.sums=append(ts.sums,sum(ts.chr10$gamma))
ts.sums=append(ts.sums,sum(ts.chr11$gamma))
ts.sums=append(ts.sums,sum(ts.chr12$gamma))
ts.sums=append(ts.sums,sum(ts.chr13$gamma))
ts.sums=append(ts.sums,sum(ts.chr14$gamma))
ts.sums=append(ts.sums,sum(ts.chr15$gamma))
ts.sums=append(ts.sums,sum(ts.chr17$gamma))
ts.sums=append(ts.sums,sum(ts.chr18$gamma))
ts.sums=append(ts.sums,sum(ts.chr19$gamma))
ts.sums=append(ts.sums,sum(ts.chr20$gamma))
ts.sums=append(ts.sums,sum(ts.chr21$gamma))
ts.sums=append(ts.sums,sum(ts.chr22$gamma))
ts.sums=append(ts.sums,sum(ts.chr23$gamma))
ts.sums=append(ts.sums,sum(ts.chr24$gamma))
ts.sums=append(ts.sums,sum(ts.chr25$gamma))
ts.sums=append(ts.sums,sum(ts.chr26$gamma))
ts.sums=append(ts.sums,sum(ts.chr27$gamma))
ts.sums=append(ts.sums,sum(ts.chr28$gamma))
ts.sums=append(ts.sums,sum(ts.chr29$gamma))
ts.sums=append(ts.sums,sum(ts.chr31$gamma))
ts.sums=append(ts.sums,sum(ts.chr32$gamma))
ts.sums=append(ts.sums,sum(ts.chr33$gamma))
ts.sums=append(ts.sums,sum(ts.chr34$gamma))
ts.sums=append(ts.sums,sum(ts.chr35$gamma))
ts.sums=append(ts.sums,sum(ts.chr36$gamma))
ts.sums=append(ts.sums,sum(ts.chr37$gamma))
ts.sums=append(ts.sums,sum(ts.chrz$gamma))


### Plot sum of PIP > 0.01 vs chromosome length-----------------------------

# Make sure to first filter by PIP when reading in the parameter data, then
# rerun the commands above!

lengths <- c(119023421,76187387,156035725,116801625,73257097,9617204,63258489,36085389,38459648,31262510,25880253,20272128,21491857,20890524,18810845,16541138,13985943,11194341,12073725,11382101,15277844,7507825,5297670,6778862,7098401,2102120,6843954,5236451,5553549,1648998,1590086,784579,437724,606149,523230,338027,276370,90132487)

par(mfrow=c(1,2))
plot(lengths,vc.sums,pch=20,xlab='Chromosome Length (Mb)',ylab='Sum PIP >= 0.01',col='grey')
abline(lm(vc.sums~lengths),lty=2)
points(76187387,vc.sums[2],pch=20,col='red')
points(90132487,vc.sums[38],pch=20,col='red')

plot(lengths,ts.sums,pch=20,xlab='Chromosome Length (Mb)',ylab='Sum PIP >= 0.01',col='grey')
abline(lm(ts.sums~lengths),lty=2)
points(156035725,ts.sums[3],pch=20,col='blue')

plot(vc.sums,pch=20)
plot(ts.sums,pch=20)

### Plot mean +/- SE PIP for each chromosome--------------------------------

## Output at 3 x 8

par(mfrow=c(1,1))
counter <- 1
plot(1,type='n',xlim=c(0,37),ylim=c(0,0.0000225),ylab='PIP',xlab='Chromosome')
for (c in vc.means){
  print(vc.maxs[counter])
  points(x=counter-1,pch=20,y=vc.means[counter],col='grey')
  segments(x0=counter-1,x1=counter-1,y0=vc.means[counter]-vc.stds[counter],y1=vc.means[counter]+vc.stds[counter],lwd=1.5)
  counter <- counter + 1
}

counter <- 1
plot(1,type='n',xlim=c(0,37),ylim=c(0,0.6),ylab='PIP',xlab='Chromosome')
for (c in vc.maxs){
  print(vc.maxs[counter])
  points(x=counter-1,pch=20,y=vc.maxs[counter],col='grey')
  counter <- counter + 1
}

counter <- 1
plot(1,type='n',xlim=c(0,37),ylim=c(0,1),ylab='Sum PIP >= 0.01',xlab='Chromosome')
for (c in vc.sums){
  print(vc.sums[counter])
  points(x=counter-1,pch=20,y=vc.sums[counter],col='grey')
  counter <- counter + 1
}

counter <- 1
plot(1,type='n',xlim=c(0,37),ylim=c(0,0.00001),ylab='PIP',xlab='Chromosome')
for (c in ts.means){
  print(ts.maxs[counter])
  points(x=counter-1,pch=20,y=ts.means[counter],col='grey')
  segments(x0=counter-1,x1=counter-1,y0=ts.means[counter]-ts.stds[counter],y1=ts.means[counter]+ts.stds[counter],lwd=1.5)
  counter <- counter + 1
}

counter <- 1
plot(1,type='n',xlim=c(0,37),ylim=c(0,0.45),ylab='PIP',xlab='Chromosome')
for (c in ts.maxs){
  print(ts.maxs[counter])
  points(x=counter-1,pch=20,y=ts.maxs[counter],col='grey')
  counter <- counter + 1
}

counter <- 1
plot(1,type='n',xlim=c(0,37),ylim=c(0,1.05),ylab='Sum PIP >= 0.01',xlab='Chromosome')
for (c in ts.sums){
  print(ts.sums[counter])
  points(x=counter-1,pch=20,y=ts.sums[counter],col='grey')
  counter <- counter + 1
}
