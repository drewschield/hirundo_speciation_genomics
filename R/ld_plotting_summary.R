############################################################################
# Barn swallow linkage disequilibrium in candidate and background regions
############################################################################

# We've calculated linkage disequilibrium (r2) between SNPs in candidate and
# background genomic regions in matched bins of parental population allele
# frequency differences (see details in README and ld_allele-freq-diff.R).

# Here, we'll compare the distributions of r2 for allele frequency
# difference-matched (AFD) candidate and background SNPs.

# This script contains commands for reading in r2 values for sets of SNPs in
# bins by parental population allele frequency differences and for summarizing
# and plotting distributions of r2.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/ld')

library(data.table)
library(tidyverse)
library(ggplot2)
library(scales)

### Read in data------------------------------------------------------------

# Hybrid zones

r2.rt.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.rustica-tytleri.interchrom.hap.ld',header=T)

r2.rt.back.2 <- read.table('./r2/r2.background.AFD-bin2.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.25 <- read.table('./r2/r2.background.AFD-bin25.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.3 <- read.table('./r2/r2.background.AFD-bin3.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.35 <- read.table('./r2/r2.background.AFD-bin35.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.4 <- read.table('./r2/r2.background.AFD-bin4.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.45 <- read.table('./r2/r2.background.AFD-bin45.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.5 <- read.table('./r2/r2.background.AFD-bin5.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.55 <- read.table('./r2/r2.background.AFD-bin55.rustica-tytleri.interchrom.hap.ld',header=T)
r2.rt.back.6 <- read.table('./r2/r2.background.AFD-bin6.rustica-tytleri.interchrom.hap.ld',header=T)

r2.rg.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.rustica-gutturalis.interchrom.hap.ld',header=T)

r2.rg.back.2 <- read.table('./r2/r2.background.AFD-bin2.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.25 <- read.table('./r2/r2.background.AFD-bin25.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.3 <- read.table('./r2/r2.background.AFD-bin3.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.35 <- read.table('./r2/r2.background.AFD-bin35.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.4 <- read.table('./r2/r2.background.AFD-bin4.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.45 <- read.table('./r2/r2.background.AFD-bin45.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.5 <- read.table('./r2/r2.background.AFD-bin5.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.55 <- read.table('./r2/r2.background.AFD-bin55.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.rg.back.6 <- read.table('./r2/r2.background.AFD-bin6.rustica-gutturalis.interchrom.hap.ld',header=T)

r2.tg.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.tytleri-gutturalis.interchrom.hap.ld',header=T)

r2.tg.back.2 <- read.table('./r2/r2.background.AFD-bin2.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.25 <- read.table('./r2/r2.background.AFD-bin25.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.3 <- read.table('./r2/r2.background.AFD-bin3.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.35 <- read.table('./r2/r2.background.AFD-bin35.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.4 <- read.table('./r2/r2.background.AFD-bin4.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.45 <- read.table('./r2/r2.background.AFD-bin45.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.5 <- read.table('./r2/r2.background.AFD-bin5.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.55 <- read.table('./r2/r2.background.AFD-bin55.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.tg.back.6 <- read.table('./r2/r2.background.AFD-bin6.tytleri-gutturalis.interchrom.hap.ld',header=T)

# Parentals

r2.ru.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.rustica.interchrom.hap.ld',header=T)
r2.ru.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.rustica.interchrom.hap.ld',header=T)

r2.ru.back.2 <- read.table('./r2/r2.background.AFD-bin2.rustica.interchrom.hap.ld',header=T)
r2.ru.back.25 <- read.table('./r2/r2.background.AFD-bin25.rustica.interchrom.hap.ld',header=T)
r2.ru.back.3 <- read.table('./r2/r2.background.AFD-bin3.rustica.interchrom.hap.ld',header=T)
r2.ru.back.35 <- read.table('./r2/r2.background.AFD-bin35.rustica.interchrom.hap.ld',header=T)
r2.ru.back.4 <- read.table('./r2/r2.background.AFD-bin4.rustica.interchrom.hap.ld',header=T)
r2.ru.back.45 <- read.table('./r2/r2.background.AFD-bin45.rustica.interchrom.hap.ld',header=T)
r2.ru.back.5 <- read.table('./r2/r2.background.AFD-bin5.rustica.interchrom.hap.ld',header=T)
r2.ru.back.55 <- read.table('./r2/r2.background.AFD-bin55.rustica.interchrom.hap.ld',header=T)
r2.ru.back.6 <- read.table('./r2/r2.background.AFD-bin6.rustica.interchrom.hap.ld',header=T)

r2.ty.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.tytleri.interchrom.hap.ld',header=T)
r2.ty.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.tytleri.interchrom.hap.ld',header=T)

r2.ty.back.2 <- read.table('./r2/r2.background.AFD-bin2.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.25 <- read.table('./r2/r2.background.AFD-bin25.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.3 <- read.table('./r2/r2.background.AFD-bin3.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.35 <- read.table('./r2/r2.background.AFD-bin35.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.4 <- read.table('./r2/r2.background.AFD-bin4.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.45 <- read.table('./r2/r2.background.AFD-bin45.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.5 <- read.table('./r2/r2.background.AFD-bin5.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.55 <- read.table('./r2/r2.background.AFD-bin55.tytleri.interchrom.hap.ld',header=T)
r2.ty.back.6 <- read.table('./r2/r2.background.AFD-bin6.tytleri.interchrom.hap.ld',header=T)

r2.gu.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.gutturalis.interchrom.hap.ld',header=T)
r2.gu.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.gutturalis.interchrom.hap.ld',header=T)

r2.gu.back.2 <- read.table('./r2/r2.background.AFD-bin2.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.25 <- read.table('./r2/r2.background.AFD-bin25.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.3 <- read.table('./r2/r2.background.AFD-bin3.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.35 <- read.table('./r2/r2.background.AFD-bin35.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.4 <- read.table('./r2/r2.background.AFD-bin4.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.45 <- read.table('./r2/r2.background.AFD-bin45.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.5 <- read.table('./r2/r2.background.AFD-bin5.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.55 <- read.table('./r2/r2.background.AFD-bin55.gutturalis.interchrom.hap.ld',header=T)
r2.gu.back.6 <- read.table('./r2/r2.background.AFD-bin6.gutturalis.interchrom.hap.ld',header=T)

# Candidate loci in simulated hybrids

r2.sim.rt.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.sim.rustica-tytleri.interchrom.hap.ld',header=T)
r2.sim.rt.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.sim.rustica-tytleri.interchrom.hap.ld',header=T)

r2.sim.rg.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.sim.rustica-gutturalis.interchrom.hap.ld',header=T)
r2.sim.rg.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.sim.rustica-gutturalis.interchrom.hap.ld',header=T)

r2.sim.tg.cand.2 <- read.table('./r2/r2.candidate.AFD-bin2.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.25 <- read.table('./r2/r2.candidate.AFD-bin25.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.3 <- read.table('./r2/r2.candidate.AFD-bin3.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.35 <- read.table('./r2/r2.candidate.AFD-bin35.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.4 <- read.table('./r2/r2.candidate.AFD-bin4.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.45 <- read.table('./r2/r2.candidate.AFD-bin45.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.5 <- read.table('./r2/r2.candidate.AFD-bin5.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.55 <- read.table('./r2/r2.candidate.AFD-bin55.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)
r2.sim.tg.cand.6 <- read.table('./r2/r2.candidate.AFD-bin6.sim.tytleri-gutturalis.interchrom.hap.ld',header=T)

### Perform comparisons of ILD at extreme AFD in trait loci / background----

wilcox.test(r2.rt.back.6$R.2,r2.rt.cand.6$R.2)
wilcox.test(r2.rg.back.6$R.2,r2.rg.cand.6$R.2)
wilcox.test(r2.tg.back.55$R.2,r2.tg.cand.55$R.2)

### Comparisons of ILD at extreme AFD for hybrids and simulated F1s---------

wilcox.test(r2.sim.rt.cand.6$R.2,r2.rt.cand.6$R.2)
wilcox.test(r2.sim.rg.cand.6$R.2,r2.rg.cand.6$R.2)
wilcox.test(r2.sim.tg.cand.55$R.2,r2.tg.cand.55$R.2)

### Calculate mean magnitude difference (effect) b/w trait/background-------

rt.mag.means <- c()
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.6$R.2,na.rm=T)/mean(r2.rt.back.6$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.55$R.2,na.rm=T)/mean(r2.rt.back.55$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.5$R.2,na.rm=T)/mean(r2.rt.back.5$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.45$R.2,na.rm=T)/mean(r2.rt.back.45$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.4$R.2,na.rm=T)/mean(r2.rt.back.4$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.35$R.2,na.rm=T)/mean(r2.rt.back.35$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.3$R.2,na.rm=T)/mean(r2.rt.back.3$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.25$R.2,na.rm=T)/mean(r2.rt.back.25$R.2,na.rm=T))
rt.mag.means = append(rt.mag.means,mean(r2.rt.cand.2$R.2,na.rm=T)/mean(r2.rt.back.2$R.2,na.rm=T))

rg.mag.means <- c()
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.6$R.2,na.rm=T)/mean(r2.rg.back.6$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.55$R.2,na.rm=T)/mean(r2.rg.back.55$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.5$R.2,na.rm=T)/mean(r2.rg.back.5$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.45$R.2,na.rm=T)/mean(r2.rg.back.45$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.4$R.2,na.rm=T)/mean(r2.rg.back.4$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.35$R.2,na.rm=T)/mean(r2.rg.back.35$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.3$R.2,na.rm=T)/mean(r2.rg.back.3$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.25$R.2,na.rm=T)/mean(r2.rg.back.25$R.2,na.rm=T))
rg.mag.means = append(rg.mag.means,mean(r2.rg.cand.2$R.2,na.rm=T)/mean(r2.rg.back.2$R.2,na.rm=T))

tg.mag.means <- c()
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.55$R.2,na.rm=T)/mean(r2.tg.back.55$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.5$R.2,na.rm=T)/mean(r2.tg.back.5$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.45$R.2,na.rm=T)/mean(r2.tg.back.45$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.4$R.2,na.rm=T)/mean(r2.tg.back.4$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.35$R.2,na.rm=T)/mean(r2.tg.back.35$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.3$R.2,na.rm=T)/mean(r2.tg.back.3$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.25$R.2,na.rm=T)/mean(r2.tg.back.25$R.2,na.rm=T))
tg.mag.means = append(tg.mag.means,mean(r2.tg.cand.2$R.2,na.rm=T)/mean(r2.tg.back.2$R.2,na.rm=T))

mean(rt.mag.means)
mean(rg.mag.means)
mean(tg.mag.means)

sd(rt.mag.means)
sd(rg.mag.means)
sd(tg.mag.means)

std <- function(x) sd(x)/sqrt(length(x))

std(rt.mag.means)
std(rg.mag.means)
std(tg.mag.means)

### Set means in increasing AFD intervals (for simple plots)----------------

rt.cand.means <- c()
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.2$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.25$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.3$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.35$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.4$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.45$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.5$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.55$R.2,na.rm=T))
rt.cand.means=append(rt.cand.means,mean(r2.rt.cand.6$R.2,na.rm=T))

rt.back.means <- c()
rt.back.means=append(rt.back.means,mean(r2.rt.back.2$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.25$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.3$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.35$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.4$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.45$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.5$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.55$R.2,na.rm=T))
rt.back.means=append(rt.back.means,mean(r2.rt.back.6$R.2,na.rm=T))

rg.cand.means <- c()
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.2$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.25$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.3$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.35$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.4$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.45$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.5$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.55$R.2,na.rm=T))
rg.cand.means=append(rg.cand.means,mean(r2.rg.cand.6$R.2,na.rm=T))

rg.back.means <- c()
rg.back.means=append(rg.back.means,mean(r2.rg.back.2$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.25$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.3$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.35$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.4$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.45$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.5$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.55$R.2,na.rm=T))
rg.back.means=append(rg.back.means,mean(r2.rg.back.6$R.2,na.rm=T))

tg.cand.means <- c()
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.2$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.25$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.3$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.35$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.4$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.45$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.5$R.2,na.rm=T))
tg.cand.means=append(tg.cand.means,mean(r2.tg.cand.55$R.2,na.rm=T))

tg.back.means <- c()
tg.back.means=append(tg.back.means,mean(r2.tg.back.2$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.25$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.3$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.35$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.4$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.45$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.5$R.2,na.rm=T))
tg.back.means=append(tg.back.means,mean(r2.tg.back.55$R.2,na.rm=T))

ru.cand.means <- c()
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.2$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.25$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.3$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.35$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.4$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.45$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.5$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.55$R.2,na.rm=T))
ru.cand.means=append(ru.cand.means,mean(r2.ru.cand.6$R.2,na.rm=T))

ru.back.means <- c()
ru.back.means=append(ru.back.means,mean(r2.ru.back.2$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.25$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.3$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.35$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.4$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.45$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.5$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.55$R.2,na.rm=T))
ru.back.means=append(ru.back.means,mean(r2.ru.back.6$R.2,na.rm=T))

ty.cand.means <- c()
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.2$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.25$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.3$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.35$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.4$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.45$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.5$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.55$R.2,na.rm=T))
ty.cand.means=append(ty.cand.means,mean(r2.ty.cand.6$R.2,na.rm=T))

ty.back.means <- c()
ty.back.means=append(ty.back.means,mean(r2.ty.back.2$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.25$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.3$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.35$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.4$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.45$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.5$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.55$R.2,na.rm=T))
ty.back.means=append(ty.back.means,mean(r2.ty.back.6$R.2,na.rm=T))

gu.cand.means <- c()
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.2$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.25$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.3$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.35$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.4$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.45$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.5$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.55$R.2,na.rm=T))
gu.cand.means=append(gu.cand.means,mean(r2.gu.cand.6$R.2,na.rm=T))

gu.back.means <- c()
gu.back.means=append(gu.back.means,mean(r2.gu.back.2$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.25$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.3$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.35$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.4$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.45$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.5$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.55$R.2,na.rm=T))
gu.back.means=append(gu.back.means,mean(r2.gu.back.6$R.2,na.rm=T))


### Plot simple splines comparing background to trait loci------------------

par(mfrow=c(2,3))
rt.cand.s <- smooth.spline(rt.cand.means,spar=0.6)
rt.back.s <- smooth.spline(rt.back.means,spar=0.6)
plot(rt.cand.s,type='l',lwd=2,col='#E28026',ylim=c(0.02,0.12),xlim=c(1,9))
lines(rt.back.s,lwd=2,col='grey')

rg.cand.s <- smooth.spline(rg.cand.means,spar=0.6)
rg.back.s <- smooth.spline(rg.back.means,spar=0.6)
plot(rg.cand.s,type='l',lwd=2,col='#8362AA',ylim=c(0.02,0.12),xlim=c(1,9))
lines(rg.back.s,lwd=2,col='grey')

tg.cand.s <- smooth.spline(tg.cand.means,spar=0.6)
tg.back.s <- smooth.spline(tg.back.means,spar=0.6)
plot(tg.cand.s,type='l',lwd=2,col='#70C6A8',ylim=c(0.02,0.12),xlim=c(1,9))
lines(tg.back.s,lwd=2,col='grey')

ru.cand.s <- smooth.spline(ru.cand.means,spar=0.6)
ru.back.s <- smooth.spline(ru.back.means,spar=0.6)
plot(ru.cand.s,type='l',lwd=2,col='#B03160',ylim=c(0.02,0.12),xlim=c(1,9))
lines(ru.back.s,lwd=2,col='grey')

ty.cand.s <- smooth.spline(ty.cand.means,spar=0.6)
ty.back.s <- smooth.spline(ty.back.means,spar=0.6)
plot(ty.cand.s,type='l',lwd=2,col='#EBB320',ylim=c(0.02,0.12),xlim=c(1,9))
lines(ty.back.s,lwd=2,col='grey')

gu.cand.s <- smooth.spline(gu.cand.means,spar=0.6)
gu.back.s <- smooth.spline(gu.back.means,spar=0.6)
plot(gu.cand.s,type='l',lwd=2,col='#6BA5CC',ylim=c(0.02,0.12),xlim=c(1,9))
lines(gu.back.s,lwd=2,col='grey')

### Plot full r2 distributions in AFD bins-----------------------------------

# Nice and clean gg boxplots
library(cowplot)

rt <- ggplot() +
  geom_boxplot(data=r2.rt.back.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.rt.back.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#E28026') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

rg <- ggplot() +
  geom_boxplot(data=r2.rg.back.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.rg.back.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#8362AA') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

tg <- ggplot() +
  geom_boxplot(data=r2.tg.back.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.tg.back.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_point(aes(x=17,y=0.2)) + 
  geom_point(aes(x=18,y=0.2)) + 
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

ru <- ggplot() +
  geom_boxplot(data=r2.ru.back.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#B03160') +
  geom_boxplot(data=r2.ru.back.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ru.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#B03160') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

ty <- ggplot() +
  geom_boxplot(data=r2.ty.back.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#EBB320') +
  geom_boxplot(data=r2.ty.back.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.ty.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#EBB320') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

gu <- ggplot() +
  geom_boxplot(data=r2.gu.back.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  geom_boxplot(data=r2.gu.back.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.gu.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#6BA5CC') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

rtsim <- ggplot() +
  geom_boxplot(data=r2.sim.rt.cand.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#E28026') +
  geom_boxplot(data=r2.sim.rt.cand.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rt.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#E28026') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

rgsim <- ggplot() +
  geom_boxplot(data=r2.sim.rg.cand.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#8362AA') +
  geom_boxplot(data=r2.sim.rg.cand.6,aes(x=17,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.rg.cand.6,aes(x=18,y=R.2),outlier.shape=NA,colour='#8362AA') +
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

tgsim <- ggplot() +
  geom_boxplot(data=r2.sim.tg.cand.2,aes(x=1,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.2,aes(x=2,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.25,aes(x=3,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.25,aes(x=4,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.3,aes(x=5,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.3,aes(x=6,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.35,aes(x=7,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.35,aes(x=8,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.4,aes(x=9,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.4,aes(x=10,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.45,aes(x=11,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.45,aes(x=12,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.5,aes(x=13,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.5,aes(x=14,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_boxplot(data=r2.sim.tg.cand.55,aes(x=15,y=R.2),outlier.shape=NA,colour='grey') +
  geom_boxplot(data=r2.tg.cand.55,aes(x=16,y=R.2),outlier.shape=NA,colour='#70C6A8') +
  geom_point(aes(x=17,y=0.2)) + 
  geom_point(aes(x=18,y=0.2)) + 
  xlab('Allele Frequency Difference') +
  ylab('Inter-chromosomal LD (r2)') +
  ylim(0,0.4) +
  theme_classic()

plot_grid(rt,rg,tg,ru,ty,gu,nrow=2,ncol=3)
plot_grid(rtsim,rgsim,tgsim,nrow=1,ncol=3)

### Plot full r2 distributions in AFD bins-----------------------------------

# Hybrid zones; background vs trait loci
boxplot(r2.rt.back.2$R.2,r2.rt.cand.2$R.2,
        r2.rt.back.25$R.2,r2.rt.cand.25$R.2,
        r2.rt.back.3$R.2,r2.rt.cand.3$R.2,
        r2.rt.back.35$R.2,r2.rt.cand.35$R.2,
        r2.rt.back.4$R.2,r2.rt.cand.4$R.2,
        r2.rt.back.45$R.2,r2.rt.cand.45$R.2,
        r2.rt.back.5$R.2,r2.rt.cand.5$R.2,
        r2.rt.back.55$R.2,r2.rt.cand.55$R.2,
        r2.rt.back.6$R.2,r2.rt.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

boxplot(r2.rg.back.2$R.2,r2.rg.cand.2$R.2,
        r2.rg.back.25$R.2,r2.rg.cand.25$R.2,
        r2.rg.back.3$R.2,r2.rg.cand.3$R.2,
        r2.rg.back.35$R.2,r2.rg.cand.35$R.2,
        r2.rg.back.4$R.2,r2.rg.cand.4$R.2,
        r2.rg.back.45$R.2,r2.rg.cand.45$R.2,
        r2.rg.back.5$R.2,r2.rg.cand.5$R.2,
        r2.rg.back.55$R.2,r2.rg.cand.55$R.2,
        r2.rg.back.6$R.2,r2.rg.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

boxplot(r2.tg.back.2$R.2,r2.tg.cand.2$R.2,
        r2.tg.back.25$R.2,r2.tg.cand.25$R.2,
        r2.tg.back.3$R.2,r2.tg.cand.3$R.2,
        r2.tg.back.35$R.2,r2.tg.cand.35$R.2,
        r2.tg.back.4$R.2,r2.tg.cand.4$R.2,
        r2.tg.back.45$R.2,r2.tg.cand.45$R.2,
        r2.tg.back.5$R.2,r2.tg.cand.5$R.2,
        r2.tg.back.55$R.2,r2.tg.cand.55$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

# Parentals; background vs trait loci
boxplot(r2.ru.back.2$R.2,r2.ru.cand.2$R.2,
        r2.ru.back.25$R.2,r2.ru.cand.25$R.2,
        r2.ru.back.3$R.2,r2.ru.cand.3$R.2,
        r2.ru.back.35$R.2,r2.ru.cand.35$R.2,
        r2.ru.back.4$R.2,r2.ru.cand.4$R.2,
        r2.ru.back.45$R.2,r2.ru.cand.45$R.2,
        r2.ru.back.5$R.2,r2.ru.cand.5$R.2,
        r2.ru.back.55$R.2,r2.ru.cand.55$R.2,
        r2.ru.back.6$R.2,r2.ru.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

boxplot(r2.ty.back.2$R.2,r2.ty.cand.2$R.2,
        r2.ty.back.25$R.2,r2.ty.cand.25$R.2,
        r2.ty.back.3$R.2,r2.ty.cand.3$R.2,
        r2.ty.back.35$R.2,r2.ty.cand.35$R.2,
        r2.ty.back.4$R.2,r2.ty.cand.4$R.2,
        r2.ty.back.45$R.2,r2.ty.cand.45$R.2,
        r2.ty.back.5$R.2,r2.ty.cand.5$R.2,
        r2.ty.back.55$R.2,r2.ty.cand.55$R.2,
        r2.ty.back.6$R.2,r2.ty.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))


boxplot(r2.gu.back.2$R.2,r2.gu.cand.2$R.2,
        r2.gu.back.25$R.2,r2.gu.cand.25$R.2,
        r2.gu.back.3$R.2,r2.gu.cand.3$R.2,
        r2.gu.back.35$R.2,r2.gu.cand.35$R.2,
        r2.gu.back.4$R.2,r2.gu.cand.4$R.2,
        r2.gu.back.45$R.2,r2.gu.cand.45$R.2,
        r2.gu.back.5$R.2,r2.gu.cand.5$R.2,
        r2.gu.back.55$R.2,r2.gu.cand.55$R.2,
        r2.gu.back.6$R.2,r2.gu.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

# Hybrid zones; traits in simulated F1s vs trait loci in real hybrids
boxplot(r2.sim.rt.cand.2$R.2,r2.rt.cand.2$R.2,
        r2.sim.rt.cand.25$R.2,r2.rt.cand.25$R.2,
        r2.sim.rt.cand.3$R.2,r2.rt.cand.3$R.2,
        r2.sim.rt.cand.35$R.2,r2.rt.cand.35$R.2,
        r2.sim.rt.cand.4$R.2,r2.rt.cand.4$R.2,
        r2.sim.rt.cand.45$R.2,r2.rt.cand.45$R.2,
        r2.sim.rt.cand.5$R.2,r2.rt.cand.5$R.2,
        r2.sim.rt.cand.55$R.2,r2.rt.cand.55$R.2,
        r2.sim.rt.cand.6$R.2,r2.rt.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

boxplot(r2.sim.rg.cand.2$R.2,r2.rg.cand.2$R.2,
        r2.sim.rg.cand.25$R.2,r2.rg.cand.25$R.2,
        r2.sim.rg.cand.3$R.2,r2.rg.cand.3$R.2,
        r2.sim.rg.cand.35$R.2,r2.rg.cand.35$R.2,
        r2.sim.rg.cand.4$R.2,r2.rg.cand.4$R.2,
        r2.sim.rg.cand.45$R.2,r2.rg.cand.45$R.2,
        r2.sim.rg.cand.5$R.2,r2.rg.cand.5$R.2,
        r2.sim.rg.cand.55$R.2,r2.rg.cand.55$R.2,
        r2.sim.rg.cand.6$R.2,r2.rg.cand.6$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

boxplot(r2.sim.tg.cand.2$R.2,r2.tg.cand.2$R.2,
        r2.sim.tg.cand.25$R.2,r2.tg.cand.25$R.2,
        r2.sim.tg.cand.3$R.2,r2.tg.cand.3$R.2,
        r2.sim.tg.cand.35$R.2,r2.tg.cand.35$R.2,
        r2.sim.tg.cand.4$R.2,r2.tg.cand.4$R.2,
        r2.sim.tg.cand.45$R.2,r2.tg.cand.45$R.2,
        r2.sim.tg.cand.5$R.2,r2.tg.cand.5$R.2,
        r2.sim.tg.cand.55$R.2,r2.tg.cand.55$R.2,
        col=c('grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine','grey','aquamarine'),
        outline=F,
        ylim=c(0,0.42))

