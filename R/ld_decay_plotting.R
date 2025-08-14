############################################################################
# Barn swallow linkage disequilibrium decay
############################################################################

# We've calculated linkage disequilibrium (r2) between physically-linked SNPs
# to examine the decay of LD as a function of physical distance.

# Here, we'll plot and summarize the decay of LD on autosomes and the Z
# chromosome for each set of parental and hybrid populations.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/ld_decay')

library(data.table)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

### Read in data-------------------------------------------------------------

decay.ru.auto <- read.table('./summary/pop.rustica.auto.hap.ld.decay.txt',header=T)
decay.ru.chrz <- read.table('./summary/pop.rustica.chrZ.hap.ld.decay.txt',header=T)
decay.ty.auto <- read.table('./summary/pop.tytleri.auto.hap.ld.decay.txt',header=T)
decay.ty.chrz <- read.table('./summary/pop.tytleri.chrZ.hap.ld.decay.txt',header=T)
decay.gu.auto <- read.table('./summary/pop.gutturalis.auto.hap.ld.decay.txt',header=T)
decay.gu.chrz <- read.table('./summary/pop.gutturalis.chrZ.hap.ld.decay.txt',header=T)

decay.rt.auto <- read.table('./summary/pop.rustica-tytleri.auto.hap.ld.decay.txt',header=T)
decay.rt.chrz <- read.table('./summary/pop.rustica-tytleri.chrZ.hap.ld.decay.txt',header=T)
decay.rg.auto <- read.table('./summary/pop.rustica-gutturalis.auto.hap.ld.decay.txt',header=T)
decay.rg.chrz <- read.table('./summary/pop.rustica-gutturalis.chrZ.hap.ld.decay.txt',header=T)
decay.tg.auto <- read.table('./summary/pop.tytleri-gutturalis.auto.hap.ld.decay.txt',header=T)
decay.tg.chrz <- read.table('./summary/pop.tytleri-gutturalis.chrZ.hap.ld.decay.txt',header=T)

### Plot LD decay------------------------------------------------------------

par(mfrow=c(2,3))
plot(decay.ru.chrz$distance,decay.ru.chrz$mean,type='l',lwd=2,ylim=c(0.02,0.2),lty=3,col='#B03160',xlab='Distance (bp)',ylab='Linkage Disequilibrium (r2)')
lines(decay.ru.auto$distance,decay.ru.auto$mean,type='l',lwd=2,col='#B03160')

plot(decay.ty.chrz$distance,decay.ty.chrz$mean,type='l',lwd=2,ylim=c(0.02,0.2),lty=3,col='#EBB320',xlab='Distance (bp)',ylab='Linkage Disequilibrium (r2)')
lines(decay.ty.auto$distance,decay.ty.auto$mean,type='l',lwd=2,col='#EBB320')

plot(decay.gu.chrz$distance,decay.gu.chrz$mean,type='l',lwd=2,ylim=c(0.02,0.2),lty=3,col='#6BA5CC',xlab='Distance (bp)',ylab='Linkage Disequilibrium (r2)')
lines(decay.gu.auto$distance,decay.gu.auto$mean,type='l',lwd=2,col='#6BA5CC')

plot(decay.rt.chrz$distance,decay.rt.chrz$mean,type='l',lwd=2,ylim=c(0.02,0.2),lty=3,col='#E28026',xlab='Distance (bp)',ylab='Linkage Disequilibrium (r2)')
lines(decay.rt.auto$distance,decay.rt.auto$mean,type='l',lwd=2,col='#E28026')

plot(decay.rg.chrz$distance,decay.rg.chrz$mean,type='l',lwd=2,ylim=c(0.02,0.2),lty=3,col='#8362AA',xlab='Distance (bp)',ylab='Linkage Disequilibrium (r2)')
lines(decay.rg.auto$distance,decay.rg.auto$mean,type='l',lwd=2,col='#8362AA')

plot(decay.tg.chrz$distance,decay.tg.chrz$mean,type='l',lwd=2,ylim=c(0.02,0.2),lty=3,col='#70C6A8',xlab='Distance (bp)',ylab='Linkage Disequilibrium (r2)')
lines(decay.tg.auto$distance,decay.tg.auto$mean,type='l',lwd=2,col='#70C6A8')
