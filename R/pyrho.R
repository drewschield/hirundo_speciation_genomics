############################################################################
# Barn swallow recombination rate variation
############################################################################

# We estimated fine-scale population recombination rate in H. rustica using
# pyrho based on a population history inferred using SMC++.

# This script contains commands for parsing, plotting, and performing summary
# statistics on genome-wide recombination rate variation.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/pyrho/')

library(data.table)
library(tidyverse)
library(scales)

options('stringsAsFactors'=FALSE)

### Read in data------------------------------------------------------------

pyrho.1mb <- read.table('./rustica.rmap.1Mb-100kb.txt',header=T,na.strings = ".")

## Parse autosomes from Z chromosome
pyrho.1mb.auto <- pyrho.1mb %>% filter(!str_detect(chrom, 'NC_053488.1'))
pyrho.1mb.chrZ <- pyrho.1mb %>% filter(str_detect(chrom, 'NC_053488.1'))

### Plot genome-wide recombination rate variation---------------------------

## Plot autosomes and Z chromosome on same x-scale (different y-axis scales)
par(mfrow=c(2,1))
plot(pyrho.1mb$rate,type='l',lwd=1.5,ylab='Recombination Rate')

place <- 0
for (scaff in unique(pyrho.1mb$chrom)){
  num <- max(as.integer(row.names(pyrho.1mb %>% filter(str_detect(scaff,chrom)))))
  place <- place + num
  #print(place)
  abline(v=place)
}

plot(pyrho.1mb.chrZ$rate,type='l',lwd=1.5,xlim=c(0,10522),ylab='Recombination Rate')

### Summary statistics------------------------------------------------------

mean(pyrho.1mb$rate,na.rm=T)
sd(pyrho.1mb$rate,na.rm=T)

mean(pyrho.1mb.auto$rate,na.rm=T)
sd(pyrho.1mb.auto$rate,na.rm=T)

mean(pyrho.1mb.chrZ$rate,na.rm=T)
sd(pyrho.1mb.chrZ$rate,na.rm=T)
