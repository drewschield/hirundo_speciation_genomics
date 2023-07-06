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
tg_ty.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/tytleri-gutturalis.locus.',args[1],'.tytleri.vcf.gz'))
tg_gu.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/tytleri-gutturalis.locus.',args[1],'.gutturalis.vcf.gz'))
tg_tg.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/tytleri-gutturalis.locus.',args[1],'.hybrids.vcf.gz'))

## Extract genotype matrices

tg_ty.bkgd.gt <- extract.gt(tg_ty.bkgd.vcfr)
tg_gu.bkgd.gt <- extract.gt(tg_gu.bkgd.vcfr)
tg_tg.bkgd.gt <- extract.gt(tg_tg.bkgd.vcfr)

## Set locus sets
tg.bkgd.loci.dat <- as.data.frame(tg_tg.bkgd.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

## Prepare data for analysis
tg_ty.bkgd.dat <- prepare.data(tg_ty.bkgd.gt, loci.data = tg.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.bkgd.gt, parental2 = tg_gu.bkgd.gt)
tg_gu.bkgd.dat <- prepare.data(tg_gu.bkgd.gt, loci.data = tg.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.bkgd.gt, parental2 = tg_gu.bkgd.gt)
tg_tg.bkgd.dat <- prepare.data(tg_tg.bkgd.gt, loci.data = tg.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.bkgd.gt, parental2 = tg_gu.bkgd.gt)

## Estimate hybrid index
tg_ty.bkgd.h <- est.h(tg_ty.bkgd.dat, loci.data = tg.bkgd.loci.dat, fixed = FALSE)
tg_gu.bkgd.h <- est.h(tg_gu.bkgd.dat, loci.data = tg.bkgd.loci.dat, fixed = FALSE)
tg_tg.bkgd.h <- est.h(tg_tg.bkgd.dat, loci.data = tg.bkgd.loci.dat, fixed = FALSE)

## Query hybrid index values and write to `data.csv`
#print(tg_ty.bkgd.h)
#print(tg_tg.bkgd.h)
#print(tg_gu.bkgd.h)

## Append data together
tg.bkgd.h <- rbind(tg_ty.bkgd.h, tg_tg.bkgd.h, tg_gu.bkgd.h)

## Write output
write.table(file=paste0('./results_introgress/tytleri-gutturalis.locus.',args[1],'.txt'), tg.bkgd.h, quote=F, row.names=F, sep="\t")

## Quit
quit(save="no")
