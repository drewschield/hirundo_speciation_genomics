############################################################################
# Barn swallow hybrid index / heterozygosity analyses (Introgress)
############################################################################

# This script contains commands for analysis of individual hybrid index and
# heterozygosity to characterize hybrid classes from the barn swallow hybrid
# zones in Asia.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/introgress/')

# Install packages
install.packages('Rmisc')
install.packages('tidyverse')
install.packages('ape')
install.packages('pegas')
install.packages('seqinr')
install.packages('ggplot2')
install.packages('adegenet')
install.packages('vcfR')
install.packages('introgress')

# Note, you may need to install introgress from source tarball: http://cran.nexr.com/web/packages/introgress/index.html

# Load libraries
library(Rmisc)
library(tidyverse)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(adegenet)
library(vcfR)
library(introgress)
library(data.table)

### Read in VCFs------------------------------------------------------------

rt_ru.vcfr <- read.vcfR('./input/anc-info.snps.rustica-tytleri_rustica.vcf.gz')
rt_ty.vcfr <- read.vcfR('./input/anc-info.snps.rustica-tytleri_tytleri.vcf.gz')
rt_rt.vcfr <- read.vcfR('./input/anc-info.snps.rustica-tytleri_hybrids.vcf.gz')

rg_ru.vcfr <- read.vcfR('./input/anc-info.snps.rustica-gutturalis_rustica.vcf.gz')
rg_gu.vcfr <- read.vcfR('./input/anc-info.snps.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.vcfr <- read.vcfR('./input/anc-info.snps.rustica-gutturalis_hybrids.vcf.gz')

tg_ty.vcfr <- read.vcfR('./input/anc-info.snps.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.vcfr <- read.vcfR('./input/anc-info.snps.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.vcfr <- read.vcfR('./input/anc-info.snps.tytleri-gutturalis_hybrids.vcf.gz')

### Extract genotype matrices-----------------------------------------------

rt_ru.gt <- extract.gt(rt_ru.vcfr)
rt_ty.gt <- extract.gt(rt_ty.vcfr)
rt_rt.gt <- extract.gt(rt_rt.vcfr)

rg_ru.gt <- extract.gt(rg_ru.vcfr)
rg_gu.gt <- extract.gt(rg_gu.vcfr)
rg_rg.gt <- extract.gt(rg_rg.vcfr)

tg_ty.gt <- extract.gt(tg_ty.vcfr)
tg_gu.gt <- extract.gt(tg_gu.vcfr)
tg_tg.gt <- extract.gt(tg_tg.vcfr)

### Set locus sets----------------------------------------------------------

rt.loci.dat <- as.data.frame(rt_rt.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.loci.dat <- as.data.frame(rg_rg.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.loci.dat <- as.data.frame(tg_tg.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

### Prepare data for analysis-----------------------------------------------

rt_ru.dat <- prepare.data(rt_ru.gt, loci.data = rt.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.gt, parental2 = rt_ty.gt)
rt_ty.dat <- prepare.data(rt_ty.gt, loci.data = rt.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.gt, parental2 = rt_ty.gt)
rt_rt.dat <- prepare.data(rt_rt.gt, loci.data = rt.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.gt, parental2 = rt_ty.gt)

rg_ru.dat <- prepare.data(rg_ru.gt, loci.data = rg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.gt, parental2 = rg_gu.gt)
rg_gu.dat <- prepare.data(rg_gu.gt, loci.data = rg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.gt, parental2 = rg_gu.gt)
rg_rg.dat <- prepare.data(rg_rg.gt, loci.data = rg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.gt, parental2 = rg_gu.gt)

tg_ty.dat <- prepare.data(tg_ty.gt, loci.data = tg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.gt, parental2 = tg_gu.gt)
tg_gu.dat <- prepare.data(tg_gu.gt, loci.data = tg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.gt, parental2 = tg_gu.gt)
tg_tg.dat <- prepare.data(tg_tg.gt, loci.data = tg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.gt, parental2 = tg_gu.gt)

### Estimate hybrid index---------------------------------------------------

rt_ru_h <- est.h(rt_ru.dat, loci.data = rt.loci.dat, fixed = FALSE)
rt_ty_h <- est.h(rt_ty.dat, loci.data = rt.loci.dat, fixed = FALSE)
rt_rt_h <- est.h(rt_rt.dat, loci.data = rt.loci.dat, fixed = FALSE)

rg_ru_h <- est.h(rg_ru.dat, loci.data = rg.loci.dat, fixed = FALSE)
rg_gu_h <- est.h(rg_gu.dat, loci.data = rg.loci.dat, fixed = FALSE)
rg_rg_h <- est.h(rg_rg.dat, loci.data = rg.loci.dat, fixed = FALSE)

tg_ty_h <- est.h(tg_ty.dat, loci.data = tg.loci.dat, fixed = FALSE)
tg_gu_h <- est.h(tg_gu.dat, loci.data = tg.loci.dat, fixed = FALSE)
tg_tg_h <- est.h(tg_tg.dat, loci.data = tg.loci.dat, fixed = FALSE)

### Calculate heterozygosity------------------------------------------------

rt_ru_het <- calc.intersp.het(rt_ru.dat)
rt_ty_het <- calc.intersp.het(rt_ty.dat)
rt_rt_het <- calc.intersp.het(rt_rt.dat)

rg_ru_het <- calc.intersp.het(rg_ru.dat)
rg_gu_het <- calc.intersp.het(rg_gu.dat)
rg_rg_het <- calc.intersp.het(rg_rg.dat)

tg_ty_het <- calc.intersp.het(tg_ty.dat)
tg_gu_het <- calc.intersp.het(tg_gu.dat)
tg_tg_het <- calc.intersp.het(tg_tg.dat)

### Plot--------------------------------------------------------------------

par(mfrow=c(1,1))
plot(rt_ru_h$h,rt_ru_het,ylim=c(0,1),xlim=c(0,1),pch=20,col='darkred')
points(rt_ty_h$h,rt_ty_het,pch=20,col='goldenrod3')
points(rt_rt_h$h,rt_rt_het,pch=20,col='darkorange')
segments(0,0, 0.5, 1, lwd=1.5, col="black")
segments(0.5,1, 1,0, lwd=1.5, col="black")

plot(rg_ru_h$h,rg_ru_het,ylim=c(0,1),xlim=c(0,1),pch=20,col='darkred')
points(rg_gu_h$h,rg_gu_het,pch=20,col='skyblue3')
points(rg_rg_h$h,rg_rg_het,pch=20,col='violet')
segments(0,0, 0.5, 1, lwd=1.5, col="black")
segments(0.5,1, 1,0, lwd=1.5, col="black")

plot(tg_ty_h$h,tg_ty_het,ylim=c(0,1),xlim=c(0,1),pch=20,col='goldenrod3')
points(tg_gu_h$h,tg_gu_het,pch=20,col='skyblue3')
points(tg_tg_h$h,tg_tg_het,pch=20,col='aquamarine3')
segments(0,0, 0.5, 1, lwd=1.5, col="black")
segments(0.5,1, 1,0, lwd=1.5, col="black")

### Output data to tables---------------------------------------------------

## Append data together
rt.h <- rbind(rt_ru_h, rt_rt_h, rt_ty_h)
rg.h <- rbind(rg_ru_h, rg_rg_h, rg_gu_h)
tg.h <- rbind(tg_ty_h, tg_tg_h, tg_gu_h)

## Write output
write.table(file='./rt_hybrid-index.txt', rt.h, quote=F, row.names=F, sep="\t")
write.table(file='./rg_hybrid-index.txt', rg.h, quote=F, row.names=F, sep="\t")
write.table(file='./tg_hybrid-index.txt', tg.h, quote=F, row.names=F, sep="\t")
