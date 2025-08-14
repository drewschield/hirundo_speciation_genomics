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
rt_ru.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/rustica-tytleri.locus.',args[1],'.rustica.vcf.gz'))
rt_ty.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/rustica-tytleri.locus.',args[1],'.tytleri.vcf.gz'))
rt_rt.bkgd.vcfr <- read.vcfR(paste0('./input_background_introgress/rustica-tytleri.locus.',args[1],'.hybrids.vcf.gz'))

## Extract genotype matrices

rt_ru.bkgd.gt <- extract.gt(rt_ru.bkgd.vcfr)
rt_ty.bkgd.gt <- extract.gt(rt_ty.bkgd.vcfr)
rt_rt.bkgd.gt <- extract.gt(rt_rt.bkgd.vcfr)

## Set locus sets
rt.bkgd.loci.dat <- as.data.frame(rt_rt.bkgd.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

## Prepare data for analysis
rt_ru.bkgd.dat <- prepare.data(rt_ru.bkgd.gt, loci.data = rt.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.bkgd.gt, parental2 = rt_ty.bkgd.gt)
rt_ty.bkgd.dat <- prepare.data(rt_ty.bkgd.gt, loci.data = rt.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.bkgd.gt, parental2 = rt_ty.bkgd.gt)
rt_rt.bkgd.dat <- prepare.data(rt_rt.bkgd.gt, loci.data = rt.bkgd.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.bkgd.gt, parental2 = rt_ty.bkgd.gt)

## Estimate hybrid index
rt_ru.bkgd.h <- est.h(rt_ru.bkgd.dat, loci.data = rt.bkgd.loci.dat, fixed = FALSE)
rt_ty.bkgd.h <- est.h(rt_ty.bkgd.dat, loci.data = rt.bkgd.loci.dat, fixed = FALSE)
rt_rt.bkgd.h <- est.h(rt_rt.bkgd.dat, loci.data = rt.bkgd.loci.dat, fixed = FALSE)

## Query hybrid index values and write to `data.csv`
#print(rt_ru.bkgd.h)
#print(rt_rt.bkgd.h)
#print(rt_ty.bkgd.h)

## Append data together
rt.bkgd.h <- rbind(rt_ru.bkgd.h, rt_rt.bkgd.h, rt_ty.bkgd.h)

## Write output
write.table(file=paste0('./results_introgress/rustica-tytleri.locus.',args[1],'.txt'), rt.bkgd.h, quote=F, row.names=F, sep="\t")

## Quit
quit(save="no")
