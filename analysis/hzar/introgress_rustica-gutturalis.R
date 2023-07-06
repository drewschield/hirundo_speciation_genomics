#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Note: this script expects a command line argument in the form of 'NC_053453.1:46296788-46306788', for example.

## Working directory and dependencies

library(tidyverse)
library(ape)
library(pegas)
library(seqinr)
library(adegenet)
library(vcfR)
library(introgress)
library(data.table)

## Read in VCF data
rg_ru.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/rustica-gutturalis.locus.',args[1],'.rustica.vcf.gz'))
rg_gu.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/rustica-gutturalis.locus.',args[1],'.gutturalis.vcf.gz'))
rg_rg.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/rustica-gutturalis.locus.',args[1],'.hybrids.vcf.gz'))

## Extract genotype matrices

rg_ru.bkgd.gt <- extract.gt(rg_ru.bkgd.vcfr)
rg_gu.bkgd.gt <- extract.gt(rg_gu.bkgd.vcfr)
rg_rg.bkgd.gt <- extract.gt(rg_rg.bkgd.vcfr)

## Set locus sets
rg.bkgd.loci.dat <- as.data.frame(rg_rg.bkgd.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

## Prepare data for analysis
rg_ru.bkgd.dat <- prepare.data(rg_ru.bkgd.gt, loci.data = rg.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.bkgd.gt, parental2 = rg_gu.bkgd.gt)
rg_gu.bkgd.dat <- prepare.data(rg_gu.bkgd.gt, loci.data = rg.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.bkgd.gt, parental2 = rg_gu.bkgd.gt)
rg_rg.bkgd.dat <- prepare.data(rg_rg.bkgd.gt, loci.data = rg.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.bkgd.gt, parental2 = rg_gu.bkgd.gt)

## Estimate hybrid index
rg_ru.bkgd.h <- est.h(rg_ru.bkgd.dat, loci.data = rg.bkgd.loci.dat, fixed = FALSE)
rg_gu.bkgd.h <- est.h(rg_gu.bkgd.dat, loci.data = rg.bkgd.loci.dat, fixed = FALSE)
rg_rg.bkgd.h <- est.h(rg_rg.bkgd.dat, loci.data = rg.bkgd.loci.dat, fixed = FALSE)

## Query hybrid index values and write to `data.csv`
#print(rg_ru.bkgd.h)
#print(rg_rg.bkgd.h)
#print(rg_gu.bkgd.h)

## Append data together
rg.bkgd.h <- rbind(rg_ru.bkgd.h, rg_rg.bkgd.h, rg_gu.bkgd.h)

## Write output
write.table(file=paste0('./results_introgress/rustica-gutturalis.locus.',args[1],'.txt'), rg.bkgd.h, quote=F, row.names=F, sep="\t")

## Quit
quit(save="no")
