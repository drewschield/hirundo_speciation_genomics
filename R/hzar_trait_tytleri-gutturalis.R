############################################################################
# Barn swallow geographic clines in hybrid zone transects
############################################################################

# This script contains commands for analysis of individual hybrid index and
# geographic clines for 'outlier' regions of the genome versus background
# to understand if genomic regions associated with mate choice traits show
# restricted gene flow (i.e., steeper/narrower clines).

# This workflow involves:
# 1. Estimation of hybrid indices from outlier regions for the tytleri-gutturalis hybrid zone
# 2. Analysis of geographic clines based on hybrid indices

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/hzar/')

#install.packages("hzar") --> had to install from source because HZAR has been removed from CRAN(?)
#install.packages('./hzar_0.2-5.tar.gz', repos = NULL, type ="source")

library(Rmisc)
library(tidyverse)
library(data.table)

library(ape)
library(pegas)
library(seqinr)
library(adegenet)
library(vcfR)
library(introgress)

library(reshape2)
library(hzar)
library(doMC)

### Load data--------------------------------------------------------------

# After running analyses below and saving HZAR data to Rdata object
load("./candidate_results_hzar/hzar_data_tytleri-gutturalis.RData")


### Introgress - estimation of hybrid index (tytleri-gutturalis)------------

## Read in VCF data

tg_ty.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.tytleri-gutturalis_hybrids.vcf.gz')

tg_ty.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.tytleri-gutturalis_tytleri.vcf.gz')
tg_gu.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.tytleri-gutturalis_gutturalis.vcf.gz')
tg_tg.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.tytleri-gutturalis_hybrids.vcf.gz')

## Extract genotype matrices

# Outlier regions
tg_ty.kitlg.gt <- extract.gt(tg_ty.kitlg.vcfr)
tg_gu.kitlg.gt <- extract.gt(tg_gu.kitlg.vcfr)
tg_tg.kitlg.gt <- extract.gt(tg_tg.kitlg.vcfr)

tg_ty.plxnc1.gt <- extract.gt(tg_ty.plxnc1.vcfr)
tg_gu.plxnc1.gt <- extract.gt(tg_gu.plxnc1.vcfr)
tg_tg.plxnc1.gt <- extract.gt(tg_tg.plxnc1.vcfr)

tg_ty.spef2.prlr.gt <- extract.gt(tg_ty.spef2.prlr.vcfr)
tg_gu.spef2.prlr.gt <- extract.gt(tg_gu.spef2.prlr.vcfr)
tg_tg.spef2.prlr.gt <- extract.gt(tg_tg.spef2.prlr.vcfr)

tg_ty.slc45a2.gt <- extract.gt(tg_ty.slc45a2.vcfr)
tg_gu.slc45a2.gt <- extract.gt(tg_gu.slc45a2.vcfr)
tg_tg.slc45a2.gt <- extract.gt(tg_tg.slc45a2.vcfr)

tg_ty.bnc2.gt <- extract.gt(tg_ty.bnc2.vcfr)
tg_gu.bnc2.gt <- extract.gt(tg_gu.bnc2.vcfr)
tg_tg.bnc2.gt <- extract.gt(tg_tg.bnc2.vcfr)

tg_ty.gnaq.gt <- extract.gt(tg_ty.gnaq.vcfr)
tg_gu.gnaq.gt <- extract.gt(tg_gu.gnaq.vcfr)
tg_tg.gnaq.gt <- extract.gt(tg_tg.gnaq.vcfr)

tg_ty.ror2.gt <- extract.gt(tg_ty.ror2.vcfr)
tg_gu.ror2.gt <- extract.gt(tg_gu.ror2.vcfr)
tg_tg.ror2.gt <- extract.gt(tg_tg.ror2.vcfr)

tg_ty.apc.camk4.gt <- extract.gt(tg_ty.apc.camk4.vcfr)
tg_gu.apc.camk4.gt <- extract.gt(tg_gu.apc.camk4.vcfr)
tg_tg.apc.camk4.gt <- extract.gt(tg_tg.apc.camk4.vcfr)

tg_ty.ice1.lncrna.gt <- extract.gt(tg_ty.ice1.lncrna.vcfr)
tg_gu.ice1.lncrna.gt <- extract.gt(tg_gu.ice1.lncrna.vcfr)
tg_tg.ice1.lncrna.gt <- extract.gt(tg_tg.ice1.lncrna.vcfr)

tg_ty.pde1c.gt <- extract.gt(tg_ty.pde1c.vcfr)
tg_gu.pde1c.gt <- extract.gt(tg_gu.pde1c.vcfr)
tg_tg.pde1c.gt <- extract.gt(tg_tg.pde1c.vcfr)

## Set locus sets

tg.kitlg.loci.dat <- as.data.frame(tg_tg.kitlg.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.plxnc1.loci.dat <- as.data.frame(tg_tg.plxnc1.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.spef2.prlr.loci.dat <- as.data.frame(tg_tg.spef2.prlr.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.slc45a2.loci.dat <- as.data.frame(tg_tg.slc45a2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.bnc2.loci.dat <- as.data.frame(tg_tg.bnc2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.gnaq.loci.dat <- as.data.frame(tg_tg.gnaq.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.ror2.loci.dat <- as.data.frame(tg_tg.ror2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.apc.camk4.loci.dat <- as.data.frame(tg_tg.apc.camk4.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.ice1.lncrna.loci.dat <- as.data.frame(tg_tg.ice1.lncrna.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

tg.pde1c.loci.dat <- as.data.frame(tg_tg.pde1c.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

## Prepare data for analysis

tg_ty.kitlg.dat <- prepare.data(tg_ty.kitlg.gt, loci.data = tg.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.kitlg.gt, parental2 = tg_gu.kitlg.gt)
tg_gu.kitlg.dat <- prepare.data(tg_gu.kitlg.gt, loci.data = tg.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.kitlg.gt, parental2 = tg_gu.kitlg.gt)
tg_tg.kitlg.dat <- prepare.data(tg_tg.kitlg.gt, loci.data = tg.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.kitlg.gt, parental2 = tg_gu.kitlg.gt)

tg_ty.plxnc1.dat <- prepare.data(tg_ty.plxnc1.gt, loci.data = tg.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.plxnc1.gt, parental2 = tg_gu.plxnc1.gt)
tg_gu.plxnc1.dat <- prepare.data(tg_gu.plxnc1.gt, loci.data = tg.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.plxnc1.gt, parental2 = tg_gu.plxnc1.gt)
tg_tg.plxnc1.dat <- prepare.data(tg_tg.plxnc1.gt, loci.data = tg.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.plxnc1.gt, parental2 = tg_gu.plxnc1.gt)

tg_ty.spef2.prlr.dat <- prepare.data(tg_ty.spef2.prlr.gt, loci.data = tg.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.spef2.prlr.gt, parental2 = tg_gu.spef2.prlr.gt)
tg_gu.spef2.prlr.dat <- prepare.data(tg_gu.spef2.prlr.gt, loci.data = tg.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.spef2.prlr.gt, parental2 = tg_gu.spef2.prlr.gt)
tg_tg.spef2.prlr.dat <- prepare.data(tg_tg.spef2.prlr.gt, loci.data = tg.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.spef2.prlr.gt, parental2 = tg_gu.spef2.prlr.gt)

tg_ty.slc45a2.dat <- prepare.data(tg_ty.slc45a2.gt, loci.data = tg.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.slc45a2.gt, parental2 = tg_gu.slc45a2.gt)
tg_gu.slc45a2.dat <- prepare.data(tg_gu.slc45a2.gt, loci.data = tg.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.slc45a2.gt, parental2 = tg_gu.slc45a2.gt)
tg_tg.slc45a2.dat <- prepare.data(tg_tg.slc45a2.gt, loci.data = tg.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.slc45a2.gt, parental2 = tg_gu.slc45a2.gt)

tg_ty.bnc2.dat <- prepare.data(tg_ty.bnc2.gt, loci.data = tg.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.bnc2.gt, parental2 = tg_gu.bnc2.gt)
tg_gu.bnc2.dat <- prepare.data(tg_gu.bnc2.gt, loci.data = tg.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.bnc2.gt, parental2 = tg_gu.bnc2.gt)
tg_tg.bnc2.dat <- prepare.data(tg_tg.bnc2.gt, loci.data = tg.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.bnc2.gt, parental2 = tg_gu.bnc2.gt)

tg_ty.gnaq.dat <- prepare.data(tg_ty.gnaq.gt, loci.data = tg.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.gnaq.gt, parental2 = tg_gu.gnaq.gt)
tg_gu.gnaq.dat <- prepare.data(tg_gu.gnaq.gt, loci.data = tg.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.gnaq.gt, parental2 = tg_gu.gnaq.gt)
tg_tg.gnaq.dat <- prepare.data(tg_tg.gnaq.gt, loci.data = tg.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.gnaq.gt, parental2 = tg_gu.gnaq.gt)

tg_ty.ror2.dat <- prepare.data(tg_ty.ror2.gt, loci.data = tg.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.ror2.gt, parental2 = tg_gu.ror2.gt)
tg_gu.ror2.dat <- prepare.data(tg_gu.ror2.gt, loci.data = tg.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.ror2.gt, parental2 = tg_gu.ror2.gt)
tg_tg.ror2.dat <- prepare.data(tg_tg.ror2.gt, loci.data = tg.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.ror2.gt, parental2 = tg_gu.ror2.gt)

tg_ty.apc.camk4.dat <- prepare.data(tg_ty.apc.camk4.gt, loci.data = tg.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.apc.camk4.gt, parental2 = tg_gu.apc.camk4.gt)
tg_gu.apc.camk4.dat <- prepare.data(tg_gu.apc.camk4.gt, loci.data = tg.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.apc.camk4.gt, parental2 = tg_gu.apc.camk4.gt)
tg_tg.apc.camk4.dat <- prepare.data(tg_tg.apc.camk4.gt, loci.data = tg.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.apc.camk4.gt, parental2 = tg_gu.apc.camk4.gt)

tg_ty.ice1.lncrna.dat <- prepare.data(tg_ty.ice1.lncrna.gt, loci.data = tg.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.ice1.lncrna.gt, parental2 = tg_gu.ice1.lncrna.gt)
tg_gu.ice1.lncrna.dat <- prepare.data(tg_gu.ice1.lncrna.gt, loci.data = tg.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.ice1.lncrna.gt, parental2 = tg_gu.ice1.lncrna.gt)
tg_tg.ice1.lncrna.dat <- prepare.data(tg_tg.ice1.lncrna.gt, loci.data = tg.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.ice1.lncrna.gt, parental2 = tg_gu.ice1.lncrna.gt)

tg_ty.pde1c.dat <- prepare.data(tg_ty.pde1c.gt, loci.data = tg.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.pde1c.gt, parental2 = tg_gu.pde1c.gt)
tg_gu.pde1c.dat <- prepare.data(tg_gu.pde1c.gt, loci.data = tg.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.pde1c.gt, parental2 = tg_gu.pde1c.gt)
tg_tg.pde1c.dat <- prepare.data(tg_tg.pde1c.gt, loci.data = tg.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = tg_ty.pde1c.gt, parental2 = tg_gu.pde1c.gt)

## Estimate hybrid index

tg_ty.kitlg.h <- est.h(tg_ty.kitlg.dat, loci.data = tg.kitlg.loci.dat, fixed = FALSE)
tg_gu.kitlg.h <- est.h(tg_gu.kitlg.dat, loci.data = tg.kitlg.loci.dat, fixed = FALSE)
tg_tg.kitlg.h <- est.h(tg_tg.kitlg.dat, loci.data = tg.kitlg.loci.dat, fixed = FALSE)

tg_ty.plxnc1.h <- est.h(tg_ty.plxnc1.dat, loci.data = tg.plxnc1.loci.dat, fixed = FALSE)
tg_gu.plxnc1.h <- est.h(tg_gu.plxnc1.dat, loci.data = tg.plxnc1.loci.dat, fixed = FALSE)
tg_tg.plxnc1.h <- est.h(tg_tg.plxnc1.dat, loci.data = tg.plxnc1.loci.dat, fixed = FALSE)

tg_ty.spef2.prlr.h <- est.h(tg_ty.spef2.prlr.dat, loci.data = tg.spef2.prlr.loci.dat, fixed = FALSE)
tg_gu.spef2.prlr.h <- est.h(tg_gu.spef2.prlr.dat, loci.data = tg.spef2.prlr.loci.dat, fixed = FALSE)
tg_tg.spef2.prlr.h <- est.h(tg_tg.spef2.prlr.dat, loci.data = tg.spef2.prlr.loci.dat, fixed = FALSE)

tg_ty.slc45a2.h <- est.h(tg_ty.slc45a2.dat, loci.data = tg.slc45a2.loci.dat, fixed = FALSE)
tg_gu.slc45a2.h <- est.h(tg_gu.slc45a2.dat, loci.data = tg.slc45a2.loci.dat, fixed = FALSE)
tg_tg.slc45a2.h <- est.h(tg_tg.slc45a2.dat, loci.data = tg.slc45a2.loci.dat, fixed = FALSE)

tg_ty.bnc2.h <- est.h(tg_ty.bnc2.dat, loci.data = tg.bnc2.loci.dat, fixed = FALSE)
tg_gu.bnc2.h <- est.h(tg_gu.bnc2.dat, loci.data = tg.bnc2.loci.dat, fixed = FALSE)
tg_tg.bnc2.h <- est.h(tg_tg.bnc2.dat, loci.data = tg.bnc2.loci.dat, fixed = FALSE)

tg_ty.gnaq.h <- est.h(tg_ty.gnaq.dat, loci.data = tg.gnaq.loci.dat, fixed = FALSE)
tg_gu.gnaq.h <- est.h(tg_gu.gnaq.dat, loci.data = tg.gnaq.loci.dat, fixed = FALSE)
tg_tg.gnaq.h <- est.h(tg_tg.gnaq.dat, loci.data = tg.gnaq.loci.dat, fixed = FALSE)

tg_ty.ror2.h <- est.h(tg_ty.ror2.dat, loci.data = tg.ror2.loci.dat, fixed = FALSE)
tg_gu.ror2.h <- est.h(tg_gu.ror2.dat, loci.data = tg.ror2.loci.dat, fixed = FALSE)
tg_tg.ror2.h <- est.h(tg_tg.ror2.dat, loci.data = tg.ror2.loci.dat, fixed = FALSE)

tg_ty.apc.camk4.h <- est.h(tg_ty.apc.camk4.dat, loci.data = tg.apc.camk4.loci.dat, fixed = FALSE)
tg_gu.apc.camk4.h <- est.h(tg_gu.apc.camk4.dat, loci.data = tg.apc.camk4.loci.dat, fixed = FALSE)
tg_tg.apc.camk4.h <- est.h(tg_tg.apc.camk4.dat, loci.data = tg.apc.camk4.loci.dat, fixed = FALSE)

tg_ty.ice1.lncrna.h <- est.h(tg_ty.ice1.lncrna.dat, loci.data = tg.ice1.lncrna.loci.dat, fixed = FALSE)
tg_gu.ice1.lncrna.h <- est.h(tg_gu.ice1.lncrna.dat, loci.data = tg.ice1.lncrna.loci.dat, fixed = FALSE)
tg_tg.ice1.lncrna.h <- est.h(tg_tg.ice1.lncrna.dat, loci.data = tg.ice1.lncrna.loci.dat, fixed = FALSE)

tg_ty.pde1c.h <- est.h(tg_ty.pde1c.dat, loci.data = tg.pde1c.loci.dat, fixed = FALSE)
tg_gu.pde1c.h <- est.h(tg_gu.pde1c.dat, loci.data = tg.pde1c.loci.dat, fixed = FALSE)
tg_tg.pde1c.h <- est.h(tg_tg.pde1c.dat, loci.data = tg.pde1c.loci.dat, fixed = FALSE)

## Query hybrid index values and write to `data.csv`

print(tg_ty.kitlg.h)
print(tg_tg.kitlg.h)
print(tg_gu.kitlg.h)

print(tg_ty.plxnc1.h)
print(tg_tg.plxnc1.h)
print(tg_gu.plxnc1.h)

print(tg_ty.spef2.prlr.h)
print(tg_tg.spef2.prlr.h)
print(tg_gu.spef2.prlr.h)

print(tg_ty.slc45a2.h)
print(tg_tg.slc45a2.h)
print(tg_gu.slc45a2.h)

print(tg_ty.bnc2.h)
print(tg_tg.bnc2.h)
print(tg_gu.bnc2.h)

print(tg_ty.gnaq.h)
print(tg_tg.gnaq.h)
print(tg_gu.gnaq.h)

print(tg_ty.ror2.h)
print(tg_tg.ror2.h)
print(tg_gu.ror2.h)

print(tg_ty.apc.camk4.h)
print(tg_tg.apc.camk4.h)
print(tg_gu.apc.camk4.h)

print(tg_ty.ice1.lncrna.h)
print(tg_tg.ice1.lncrna.h)
print(tg_gu.ice1.lncrna.h)

print(tg_ty.pde1c.h)
print(tg_tg.pde1c.h)
print(tg_gu.pde1c.h)

### Read in and format data for HZAR (tytleri-gutturalis)--------------------

# These data include the sample, population, lat, long, locality, hybrid index
# (background), and hybrid index (outlier regions).

data.tg <- read.csv('./data.tytleri-gutturalis.csv',header=T)
names(data.tg)

## Standardize phenotype vector ranges between 0 and 1 and add to data frame
mean_rwl_norm = (data.tg$mean_rwl-min(data.tg$mean_rwl,na.rm = T))/(max(data.tg$mean_rwl,na.rm = T)-min(data.tg$mean_rwl,na.rm = T))
mean_ts_random_norm = (data.tg$mean_ts_random-min(data.tg$mean_ts_random,na.rm = T))/(max(data.tg$mean_ts_random,na.rm = T)-min(data.tg$mean_ts_random,na.rm = T))
breast_avg_bright_norm = (data.tg$breast_avg_bright-min(data.tg$breast_avg_bright,na.rm = T))/(max(data.tg$breast_avg_bright,na.rm = T)-min(data.tg$breast_avg_bright,na.rm = T))

data.tg$wing <- mean_rwl_norm
data.tg$tail <- mean_ts_random_norm
data.tg$vent <- breast_avg_bright_norm

## Subset to get locality, long, lat, HI-background, and HI-outlier
data.tg <- dplyr::select(data.tg, loc,long, lat, hi_kitlg, hi_plxnc1, hi_spef2_prlr, hi_slc45a2, hi_bnc2, hi_gnaq, hi_ror2, hi_apc_camk4, hi_ice1_lncrna, hi_pde1c)
head(data.tg)

## Calculate means, variances, and counts

data.tg.melt<-reshape2::melt(data.tg, id.vars=c("loc"))
str(data.tg.melt)
head(data.tg.melt)
data.tg.hzar<- data.tg.melt %>%
  group_by(loc, variable) %>%
  dplyr::summarise(mean=mean(value, na.rm=TRUE), var=var(value, na.rm=TRUE), count=sum(!is.na(value)))
head(data.tg.hzar)

## Change to wide format

d<-reshape2::melt(data.tg.hzar, id.vars=c("loc", "variable"))
head(d)
colnames(d)[3]<-"metric"
data.tg.hzar.cast<-dcast(d, loc~ variable + metric )

print(data.tg.hzar.cast)

## Calculate distances (km) for transect

tg.dists<-hzar.map.latLongSites(data.tg.hzar.cast$loc, data.tg.hzar.cast$lat_mean, data.tg.hzar.cast$long_mean, degrees = TRUE)

join_func<-function(y){
  y<-left_join(y, data.tg.hzar.cast, by=c("site"="loc"))
  y<-y[,-c(2:3)]
  y<-arrange(y, km)
}

## Subset points to the sampling transects

tg.dists<-filter(tg.dists, site=="zakaltoose" | site=="malamolevo" | site=="kytyleek" | site=="nikolaevska" | site=="tataurova" | site=="narin_talacha" | site=="mixed_barns" | site=="cincer_mandel_so" | site=="bulgan_soum_dorn" | site=="dashbalbar" | site=="bayan_uul" | site=="norovlin" | site=="chingis" | site=="tsenhermandal" | site=="erdene" | site=="zuunmod" | site=="beijing" | site=="shuang" | site=="qiqihar" | site=="qinhuangdao" | site=="shenyang" | site=="changchun" | site=="harbin") # keep locations in transect
print(tg.dists)

##Add km transect distances

tg.dists$km <- hzar.map.distanceFromSite(tg.dists,"kytyleek",units="Km")
tg.dists<-join_func(tg.dists)
tg.dists<-na.omit(tg.dists)
print(tg.dists)

### Set up MCMC chain parameters--------------------------------------------

## A typical chain length.  This value is the default setting in the package.
chainLength=1e5;                       

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528))


if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
  registerDoMC()
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}


### Analysis on KITLG region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$kitlg <- list();
## Space to hold the observed data
tg$kitlg$obs <- list();
## Space to hold the models to fit
tg$kitlg$models <- list();
## Space to hold the compiled fit requests
tg$kitlg$fitRs <- list();
## Space to hold the output data chains
tg$kitlg$runs <- list();
## Space to hold the analysed data
tg$kitlg$analysis <- list();

## Assign data for kitlg
tg$kitlg$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_kitlg_mean,
                             tg.dists$hi_kitlg_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$kitlg$obs);

## Make a helper function to set cline models
tg.loadkitlgmodel <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  tg$kitlg$models[[id]] <<- hzar.makeCline1DFreq(tg$kitlg$obs, scaling, tails)

tg.loadkitlgmodel("fixed","none","modelI");
tg.loadkitlgmodel("free" ,"none","modelII");
tg.loadkitlgmodel("free" ,"both","modelIII");
tg.loadkitlgmodel("free" ,"right","modelIV");
tg.loadkitlgmodel("free" ,"left","modelV");
tg.loadkitlgmodel("free" ,"mirror","modelVI");

## Check the default settings
print(tg$kitlg$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$kitlg$models <- sapply(tg$kitlg$models,
                         hzar.model.addBoxReq,
                         -30 , 2185,
                         simplify=FALSE)

## Check the updated settings
print(tg$kitlg$models)

## Compile each of the models to prepare for fitting
tg$kitlg$fitRs$init <- sapply(tg$kitlg$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=tg$kitlg$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
tg$kitlg$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$kitlg$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$kitlg$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$kitlg$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$kitlg$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$kitlg$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$kitlg$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$kitlg$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$kitlg$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$kitlg$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$kitlg$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$kitlg$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$kitlg$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$kitlg$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$kitlg$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$kitlg$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$kitlg$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$kitlg$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$kitlg$fitRs$init)

## Run the models for an initial chain
tg$kitlg$runs$init <- list()

tg$kitlg$runs$init$modelI <-
  hzar.doFit(tg$kitlg$fitRs$init$modelI)

tg$kitlg$runs$init$modelII <-
  hzar.doFit(tg$kitlg$fitRs$init$modelII)

tg$kitlg$runs$init$modelIII <-
  hzar.doFit(tg$kitlg$fitRs$init$modelIII)

tg$kitlg$runs$init$modelIV <-
  hzar.doFit(tg$kitlg$fitRs$init$modelIV)

tg$kitlg$runs$init$modelV <-
  hzar.doFit(tg$kitlg$fitRs$init$modelV)

tg$kitlg$runs$init$modelVI <-
  hzar.doFit(tg$kitlg$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$kitlg$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$kitlg$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$kitlg$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$kitlg$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$kitlg$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$kitlg$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$kitlg$fitRs$chains <-
  lapply(tg$kitlg$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$kitlg$fitRs$chains <-
  hzar.multiFitRequest(tg$kitlg$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$kitlg$runs$chains <-  hzar.doChain.multi(tg$kitlg$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$kitlg$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$kitlg$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$kitlg$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$kitlg$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$kitlg$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$kitlg$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$kitlg$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$kitlg$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$kitlg$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$kitlg$runs$init$modelI)
tg$kitlg$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$kitlg$runs$init$modelII)
tg$kitlg$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$kitlg$runs$init$modelIII)
tg$kitlg$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$kitlg$runs$init$modelIV)
tg$kitlg$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$kitlg$runs$init$modelV)
tg$kitlg$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$kitlg$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$kitlg$analysis$oDG <-
  hzar.make.obsDataGroup(tg$kitlg$analysis$initDGs)
tg$kitlg$analysis$oDG <-
  hzar.copyModelLabels(tg$kitlg$analysis$initDGs,
                       tg$kitlg$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$kitlg$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$kitlg$runs$chains,
                                hzar.dataGroup.add),
                         tg$kitlg$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$kitlg$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$kitlg$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$kitlg$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$kitlg$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$kitlg$analysis$model.name <-
        rownames(tg$kitlg$analysis$AICcTable
        )[[ which.min(tg$kitlg$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$kitlg$analysis$model.selected <-
  tg$kitlg$analysis$oDG$data.groups[[tg$kitlg$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$kitlg$analysis$model.selected,
                         names(tg$kitlg$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$kitlg$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$kitlg$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$kitlg$analysis$model.selected);
rect(xleft=tg$kitlg$analysis$model.selected$ML.cline$param.all$center-((tg$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=tg$kitlg$analysis$model.selected$ML.cline$param.all$center+((tg$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#70C6A8')

## End Analysis

#dev.off()

### Analysis on PLXNC1 region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$plxnc1 <- list();
## Space to hold the observed data
tg$plxnc1$obs <- list();
## Space to hold the models to fit
tg$plxnc1$models <- list();
## Space to hold the compiled fit requests
tg$plxnc1$fitRs <- list();
## Space to hold the output data chains
tg$plxnc1$runs <- list();
## Space to hold the analysed data
tg$plxnc1$analysis <- list();

## Assign data for plxnc1
tg$plxnc1$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_plxnc1_mean,
                             tg.dists$hi_plxnc1_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$plxnc1$obs);

## Make a helper function to set cline models
tg.loadplxnc1model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$plxnc1$models[[id]] <<- hzar.makeCline1DFreq(tg$plxnc1$obs, scaling, tails)

tg.loadplxnc1model("fixed","none","modelI");
tg.loadplxnc1model("free" ,"none","modelII");
tg.loadplxnc1model("free" ,"both","modelIII");
tg.loadplxnc1model("free" ,"right","modelIV");
tg.loadplxnc1model("free" ,"left","modelV");
tg.loadplxnc1model("free" ,"mirror","modelVI");

## Check the default settings
print(tg$plxnc1$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$plxnc1$models <- sapply(tg$plxnc1$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$plxnc1$models)

## Compile each of the models to prepare for fitting
tg$plxnc1$fitRs$init <- sapply(tg$plxnc1$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$plxnc1$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$plxnc1$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$plxnc1$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$plxnc1$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$plxnc1$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$plxnc1$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$plxnc1$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$plxnc1$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$plxnc1$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$plxnc1$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$plxnc1$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$plxnc1$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$plxnc1$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$plxnc1$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$plxnc1$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$plxnc1$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$plxnc1$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$plxnc1$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$plxnc1$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$plxnc1$fitRs$init)

## Run the models for an initial chain
tg$plxnc1$runs$init <- list()

tg$plxnc1$runs$init$modelI <-
  hzar.doFit(tg$plxnc1$fitRs$init$modelI)

tg$plxnc1$runs$init$modelII <-
  hzar.doFit(tg$plxnc1$fitRs$init$modelII)

tg$plxnc1$runs$init$modelIII <-
  hzar.doFit(tg$plxnc1$fitRs$init$modelIII)

tg$plxnc1$runs$init$modelIV <-
  hzar.doFit(tg$plxnc1$fitRs$init$modelIV)

tg$plxnc1$runs$init$modelV <-
  hzar.doFit(tg$plxnc1$fitRs$init$modelV)

tg$plxnc1$runs$init$modelVI <-
  hzar.doFit(tg$plxnc1$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$plxnc1$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$plxnc1$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$plxnc1$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$plxnc1$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$plxnc1$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$plxnc1$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$plxnc1$fitRs$chains <-
  lapply(tg$plxnc1$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$plxnc1$fitRs$chains <-
  hzar.multiFitRequest(tg$plxnc1$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$plxnc1$runs$chains <-  hzar.doChain.multi(tg$plxnc1$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$plxnc1$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$plxnc1$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$plxnc1$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$plxnc1$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$plxnc1$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$plxnc1$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$plxnc1$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$plxnc1$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$plxnc1$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$plxnc1$runs$init$modelI)
tg$plxnc1$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$plxnc1$runs$init$modelII)
tg$plxnc1$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$plxnc1$runs$init$modelIII)
tg$plxnc1$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$plxnc1$runs$init$modelIV)
tg$plxnc1$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$plxnc1$runs$init$modelV)
tg$plxnc1$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$plxnc1$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$plxnc1$analysis$oDG <-
  hzar.make.obsDataGroup(tg$plxnc1$analysis$initDGs)
tg$plxnc1$analysis$oDG <-
  hzar.copyModelLabels(tg$plxnc1$analysis$initDGs,
                       tg$plxnc1$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$plxnc1$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$plxnc1$runs$chains,
                                hzar.dataGroup.add),
                         tg$plxnc1$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$plxnc1$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$plxnc1$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$plxnc1$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$plxnc1$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$plxnc1$analysis$model.name <-
        rownames(tg$plxnc1$analysis$AICcTable
        )[[ which.min(tg$plxnc1$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$plxnc1$analysis$model.selected <-
  tg$plxnc1$analysis$oDG$data.groups[[tg$plxnc1$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$plxnc1$analysis$model.selected,
                         names(tg$plxnc1$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$plxnc1$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$plxnc1$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$plxnc1$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on SPEF2-PRLR region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$spef2.prlr <- list();
## Space to hold the observed data
tg$spef2.prlr$obs <- list();
## Space to hold the models to fit
tg$spef2.prlr$models <- list();
## Space to hold the compiled fit requests
tg$spef2.prlr$fitRs <- list();
## Space to hold the output data chains
tg$spef2.prlr$runs <- list();
## Space to hold the analysed data
tg$spef2.prlr$analysis <- list();

## Assign data for spef2.prlr
tg$spef2.prlr$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_spef2_prlr_mean,
                             tg.dists$hi_spef2_prlr_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$spef2.prlr$obs);

## Make a helper function to set cline models
tg.loadspef2.prlrmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$spef2.prlr$models[[id]] <<- hzar.makeCline1DFreq(tg$spef2.prlr$obs, scaling, tails)

tg.loadspef2.prlrmodel("fixed","none","modelI");
tg.loadspef2.prlrmodel("free" ,"none","modelII");
tg.loadspef2.prlrmodel("free" ,"both","modelIII");
tg.loadspef2.prlrmodel("free" ,"right","modelIV");
tg.loadspef2.prlrmodel("free" ,"left","modelV");
tg.loadspef2.prlrmodel("free" ,"mirror","modelVI");

## Check the default settings
print(tg$spef2.prlr$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$spef2.prlr$models <- sapply(tg$spef2.prlr$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$spef2.prlr$models)

## Compile each of the models to prepare for fitting
tg$spef2.prlr$fitRs$init <- sapply(tg$spef2.prlr$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$spef2.prlr$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$spef2.prlr$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$spef2.prlr$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$spef2.prlr$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$spef2.prlr$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$spef2.prlr$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$spef2.prlr$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$spef2.prlr$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$spef2.prlr$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$spef2.prlr$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$spef2.prlr$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$spef2.prlr$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$spef2.prlr$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$spef2.prlr$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$spef2.prlr$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$spef2.prlr$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$spef2.prlr$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$spef2.prlr$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$spef2.prlr$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$spef2.prlr$fitRs$init)

## Run the models for an initial chain
tg$spef2.prlr$runs$init <- list()

tg$spef2.prlr$runs$init$modelI <-
  hzar.doFit(tg$spef2.prlr$fitRs$init$modelI)

tg$spef2.prlr$runs$init$modelII <-
  hzar.doFit(tg$spef2.prlr$fitRs$init$modelII)

tg$spef2.prlr$runs$init$modelIII <-
  hzar.doFit(tg$spef2.prlr$fitRs$init$modelIII)

tg$spef2.prlr$runs$init$modelIV <-
  hzar.doFit(tg$spef2.prlr$fitRs$init$modelIV)

tg$spef2.prlr$runs$init$modelV <-
  hzar.doFit(tg$spef2.prlr$fitRs$init$modelV)

tg$spef2.prlr$runs$init$modelVI <-
  hzar.doFit(tg$spef2.prlr$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$spef2.prlr$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$spef2.prlr$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$spef2.prlr$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$spef2.prlr$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$spef2.prlr$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$spef2.prlr$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$spef2.prlr$fitRs$chains <-
  lapply(tg$spef2.prlr$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$spef2.prlr$fitRs$chains <-
  hzar.multiFitRequest(tg$spef2.prlr$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$spef2.prlr$runs$chains <-  hzar.doChain.multi(tg$spef2.prlr$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$spef2.prlr$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$spef2.prlr$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$spef2.prlr$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$spef2.prlr$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$spef2.prlr$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$spef2.prlr$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$spef2.prlr$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$spef2.prlr$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$spef2.prlr$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$spef2.prlr$runs$init$modelI)
tg$spef2.prlr$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$spef2.prlr$runs$init$modelII)
tg$spef2.prlr$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$spef2.prlr$runs$init$modelIII)
tg$spef2.prlr$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$spef2.prlr$runs$init$modelIV)
tg$spef2.prlr$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$spef2.prlr$runs$init$modelV)
tg$spef2.prlr$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$spef2.prlr$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$spef2.prlr$analysis$oDG <-
  hzar.make.obsDataGroup(tg$spef2.prlr$analysis$initDGs)
tg$spef2.prlr$analysis$oDG <-
  hzar.copyModelLabels(tg$spef2.prlr$analysis$initDGs,
                       tg$spef2.prlr$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$spef2.prlr$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$spef2.prlr$runs$chains,
                                hzar.dataGroup.add),
                         tg$spef2.prlr$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$spef2.prlr$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$spef2.prlr$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$spef2.prlr$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$spef2.prlr$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$spef2.prlr$analysis$model.name <-
        rownames(tg$spef2.prlr$analysis$AICcTable
        )[[ which.min(tg$spef2.prlr$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$spef2.prlr$analysis$model.selected <-
  tg$spef2.prlr$analysis$oDG$data.groups[[tg$spef2.prlr$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$spef2.prlr$analysis$model.selected,
                         names(tg$spef2.prlr$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$spef2.prlr$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$spef2.prlr$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$spef2.prlr$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on SLC45A2 region (rustica-tytleri)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$slc45a2 <- list();
## Space to hold the observed data
tg$slc45a2$obs <- list();
## Space to hold the models to fit
tg$slc45a2$models <- list();
## Space to hold the compiled fit requests
tg$slc45a2$fitRs <- list();
## Space to hold the output data chains
tg$slc45a2$runs <- list();
## Space to hold the analysed data
tg$slc45a2$analysis <- list();

## Assign data for slc45a2
tg$slc45a2$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_slc45a2_mean,
                             tg.dists$hi_slc45a2_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$slc45a2$obs);

## Make a helper function to set cline models
tg.loadslc45a2model <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  tg$slc45a2$models[[id]] <<- hzar.makeCline1DFreq(tg$slc45a2$obs, scaling, tails)

tg.loadslc45a2model("fixed","none","modelI");
tg.loadslc45a2model("free" ,"none","modelII");
tg.loadslc45a2model("free" ,"both","modelIII");
tg.loadslc45a2model("free" ,"right","modelIV");
tg.loadslc45a2model("free" ,"left","modelV");
tg.loadslc45a2model("free" ,"mirror","modelVI");

## Check the default settings
print(tg$slc45a2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$slc45a2$models <- sapply(tg$slc45a2$models,
                         hzar.model.addBoxReq,
                         -30 , 2185,
                         simplify=FALSE)

## Check the updated settings
print(tg$slc45a2$models)

## Compile each of the models to prepare for fitting
tg$slc45a2$fitRs$init <- sapply(tg$slc45a2$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=tg$slc45a2$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
tg$slc45a2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$slc45a2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$slc45a2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$slc45a2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$slc45a2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$slc45a2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$slc45a2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$slc45a2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$slc45a2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$slc45a2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$slc45a2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$slc45a2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$slc45a2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$slc45a2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$slc45a2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$slc45a2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$slc45a2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$slc45a2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$slc45a2$fitRs$init)

## Run the models for an initial chain
tg$slc45a2$runs$init <- list()

tg$slc45a2$runs$init$modelI <-
  hzar.doFit(tg$slc45a2$fitRs$init$modelI)

tg$slc45a2$runs$init$modelII <-
  hzar.doFit(tg$slc45a2$fitRs$init$modelII)

tg$slc45a2$runs$init$modelIII <-
  hzar.doFit(tg$slc45a2$fitRs$init$modelIII)

tg$slc45a2$runs$init$modelIV <-
  hzar.doFit(tg$slc45a2$fitRs$init$modelIV)

tg$slc45a2$runs$init$modelV <-
  hzar.doFit(tg$slc45a2$fitRs$init$modelV)

tg$slc45a2$runs$init$modelVI <-
  hzar.doFit(tg$slc45a2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$slc45a2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$slc45a2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$slc45a2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$slc45a2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$slc45a2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$slc45a2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$slc45a2$fitRs$chains <-
  lapply(tg$slc45a2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$slc45a2$fitRs$chains <-
  hzar.multiFitRequest(tg$slc45a2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$slc45a2$runs$chains <-  hzar.doChain.multi(tg$slc45a2$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$slc45a2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$slc45a2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$slc45a2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$slc45a2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$slc45a2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$slc45a2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$slc45a2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$slc45a2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$slc45a2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$slc45a2$runs$init$modelI)
tg$slc45a2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$slc45a2$runs$init$modelII)
tg$slc45a2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$slc45a2$runs$init$modelIII)
tg$slc45a2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$slc45a2$runs$init$modelIV)
tg$slc45a2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$slc45a2$runs$init$modelV)
tg$slc45a2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$slc45a2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$slc45a2$analysis$oDG <-
  hzar.make.obsDataGroup(tg$slc45a2$analysis$initDGs)
tg$slc45a2$analysis$oDG <-
  hzar.copyModelLabels(tg$slc45a2$analysis$initDGs,
                       tg$slc45a2$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$slc45a2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$slc45a2$runs$chains,
                                hzar.dataGroup.add),
                         tg$slc45a2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$slc45a2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$slc45a2$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$slc45a2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$slc45a2$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$slc45a2$analysis$model.name <-
        rownames(tg$slc45a2$analysis$AICcTable
        )[[ which.min(tg$slc45a2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$slc45a2$analysis$model.selected <-
  tg$slc45a2$analysis$oDG$data.groups[[tg$slc45a2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$slc45a2$analysis$model.selected,
                         names(tg$slc45a2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$slc45a2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$slc45a2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$slc45a2$analysis$model.selected);
rect(xleft=tg$slc45a2$analysis$model.selected$ML.cline$param.all$center-((tg$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=tg$slc45a2$analysis$model.selected$ML.cline$param.all$center+((tg$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#70C6A8')

## End Analysis

#dev.off()

### Analysis on BNC2 region (rustica-tytleri)-------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$bnc2 <- list();
## Space to hold the observed data
tg$bnc2$obs <- list();
## Space to hold the models to fit
tg$bnc2$models <- list();
## Space to hold the compiled fit requests
tg$bnc2$fitRs <- list();
## Space to hold the output data chains
tg$bnc2$runs <- list();
## Space to hold the analysed data
tg$bnc2$analysis <- list();

## Assign data for bnc2
tg$bnc2$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_bnc2_mean,
                             tg.dists$hi_bnc2_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$bnc2$obs);

## Make a helper function to set cline models
tg.loadbnc2model <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  tg$bnc2$models[[id]] <<- hzar.makeCline1DFreq(tg$bnc2$obs, scaling, tails)

tg.loadbnc2model("fixed","none","modelI");
tg.loadbnc2model("free" ,"none","modelII");
tg.loadbnc2model("free" ,"both","modelIII");
tg.loadbnc2model("free" ,"right","modelIV");
tg.loadbnc2model("free" ,"left","modelV");
tg.loadbnc2model("free" ,"mirror","modelVI");

## Check the default settings
print(tg$bnc2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$bnc2$models <- sapply(tg$bnc2$models,
                         hzar.model.addBoxReq,
                         -30 , 2185,
                         simplify=FALSE)

## Check the updated settings
print(tg$bnc2$models)

## Compile each of the models to prepare for fitting
tg$bnc2$fitRs$init <- sapply(tg$bnc2$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=tg$bnc2$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
tg$bnc2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$bnc2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$bnc2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$bnc2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$bnc2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$bnc2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$bnc2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$bnc2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$bnc2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$bnc2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$bnc2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$bnc2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$bnc2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$bnc2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$bnc2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$bnc2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$bnc2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$bnc2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$bnc2$fitRs$init)

## Run the models for an initial chain
tg$bnc2$runs$init <- list()

tg$bnc2$runs$init$modelI <-
  hzar.doFit(tg$bnc2$fitRs$init$modelI)

tg$bnc2$runs$init$modelII <-
  hzar.doFit(tg$bnc2$fitRs$init$modelII)

tg$bnc2$runs$init$modelIII <-
  hzar.doFit(tg$bnc2$fitRs$init$modelIII)

tg$bnc2$runs$init$modelIV <-
  hzar.doFit(tg$bnc2$fitRs$init$modelIV)

tg$bnc2$runs$init$modelV <-
  hzar.doFit(tg$bnc2$fitRs$init$modelV)

tg$bnc2$runs$init$modelVI <-
  hzar.doFit(tg$bnc2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$bnc2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$bnc2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$bnc2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$bnc2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$bnc2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$bnc2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$bnc2$fitRs$chains <-
  lapply(tg$bnc2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$bnc2$fitRs$chains <-
  hzar.multiFitRequest(tg$bnc2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$bnc2$runs$chains <-  hzar.doChain.multi(tg$bnc2$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$bnc2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$bnc2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$bnc2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$bnc2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$bnc2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$bnc2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$bnc2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$bnc2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$bnc2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$bnc2$runs$init$modelI)
tg$bnc2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$bnc2$runs$init$modelII)
tg$bnc2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$bnc2$runs$init$modelIII)
tg$bnc2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$bnc2$runs$init$modelIV)
tg$bnc2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$bnc2$runs$init$modelV)
tg$bnc2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$bnc2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$bnc2$analysis$oDG <-
  hzar.make.obsDataGroup(tg$bnc2$analysis$initDGs)
tg$bnc2$analysis$oDG <-
  hzar.copyModelLabels(tg$bnc2$analysis$initDGs,
                       tg$bnc2$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$bnc2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$bnc2$runs$chains,
                                hzar.dataGroup.add),
                         tg$bnc2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$bnc2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$bnc2$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$bnc2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$bnc2$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$bnc2$analysis$model.name <-
        rownames(tg$bnc2$analysis$AICcTable
        )[[ which.min(tg$bnc2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$bnc2$analysis$model.selected <-
  tg$bnc2$analysis$oDG$data.groups[[tg$bnc2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$bnc2$analysis$model.selected,
                         names(tg$bnc2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$bnc2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$bnc2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$bnc2$analysis$model.selected);
rect(xleft=tg$bnc2$analysis$model.selected$ML.cline$param.all$center-((tg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=tg$bnc2$analysis$model.selected$ML.cline$param.all$center+((tg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#70C6A8')

## End Analysis

#dev.off()

### Analysis on GNAQ region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$gnaq <- list();
## Space to hold the observed data
tg$gnaq$obs <- list();
## Space to hold the models to fit
tg$gnaq$models <- list();
## Space to hold the compiled fit requests
tg$gnaq$fitRs <- list();
## Space to hold the output data chains
tg$gnaq$runs <- list();
## Space to hold the analysed data
tg$gnaq$analysis <- list();

## Assign data for gnaq
tg$gnaq$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_gnaq_mean,
                             tg.dists$hi_gnaq_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$gnaq$obs);

## Make a helper function to set cline models
tg.loadgnaqmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$gnaq$models[[id]] <<- hzar.makeCline1DFreq(tg$gnaq$obs, scaling, tails)

tg.loadgnaqmodel("fixed","none","modelI");
tg.loadgnaqmodel("free" ,"none","modelII");
tg.loadgnaqmodel("free" ,"both","modelIII");
tg.loadgnaqmodel("free" ,"right","modelIV");
tg.loadgnaqmodel("free" ,"left","modelV");
tg.loadgnaqmodel("free" ,"mirror","modelVI");

## Check the default settings
print(tg$gnaq$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$gnaq$models <- sapply(tg$gnaq$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$gnaq$models)

## Compile each of the models to prepare for fitting
tg$gnaq$fitRs$init <- sapply(tg$gnaq$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$gnaq$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$gnaq$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$gnaq$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$gnaq$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$gnaq$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$gnaq$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$gnaq$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$gnaq$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$gnaq$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$gnaq$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$gnaq$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$gnaq$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$gnaq$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$gnaq$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$gnaq$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$gnaq$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$gnaq$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$gnaq$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$gnaq$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$gnaq$fitRs$init)

## Run the models for an initial chain
tg$gnaq$runs$init <- list()

tg$gnaq$runs$init$modelI <-
  hzar.doFit(tg$gnaq$fitRs$init$modelI)

tg$gnaq$runs$init$modelII <-
  hzar.doFit(tg$gnaq$fitRs$init$modelII)

tg$gnaq$runs$init$modelIII <-
  hzar.doFit(tg$gnaq$fitRs$init$modelIII)

tg$gnaq$runs$init$modelIV <-
  hzar.doFit(tg$gnaq$fitRs$init$modelIV)

tg$gnaq$runs$init$modelV <-
  hzar.doFit(tg$gnaq$fitRs$init$modelV)

tg$gnaq$runs$init$modelVI <-
  hzar.doFit(tg$gnaq$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$gnaq$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$gnaq$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$gnaq$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$gnaq$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$gnaq$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$gnaq$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$gnaq$fitRs$chains <-
  lapply(tg$gnaq$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$gnaq$fitRs$chains <-
  hzar.multiFitRequest(tg$gnaq$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$gnaq$runs$chains <-  hzar.doChain.multi(tg$gnaq$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$gnaq$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$gnaq$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$gnaq$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$gnaq$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$gnaq$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$gnaq$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$gnaq$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$gnaq$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$gnaq$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$gnaq$runs$init$modelI)
tg$gnaq$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$gnaq$runs$init$modelII)
tg$gnaq$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$gnaq$runs$init$modelIII)
tg$gnaq$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$gnaq$runs$init$modelIV)
tg$gnaq$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$gnaq$runs$init$modelV)
tg$gnaq$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$gnaq$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$gnaq$analysis$oDG <-
  hzar.make.obsDataGroup(tg$gnaq$analysis$initDGs)
tg$gnaq$analysis$oDG <-
  hzar.copyModelLabels(tg$gnaq$analysis$initDGs,
                       tg$gnaq$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$gnaq$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$gnaq$runs$chains,
                                hzar.dataGroup.add),
                         tg$gnaq$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$gnaq$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$gnaq$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$gnaq$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$gnaq$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$gnaq$analysis$model.name <-
        rownames(tg$gnaq$analysis$AICcTable
        )[[ which.min(tg$gnaq$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$gnaq$analysis$model.selected <-
  tg$gnaq$analysis$oDG$data.groups[[tg$gnaq$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$gnaq$analysis$model.selected,
                         names(tg$gnaq$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$gnaq$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$gnaq$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$gnaq$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on ROR2 region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$ror2 <- list();
## Space to hold the observed data
tg$ror2$obs <- list();
## Space to hold the models to fit
tg$ror2$models <- list();
## Space to hold the compiled fit requests
tg$ror2$fitRs <- list();
## Space to hold the output data chains
tg$ror2$runs <- list();
## Space to hold the analysed data
tg$ror2$analysis <- list();

## Assign data for ror2
tg$ror2$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_ror2_mean,
                             tg.dists$hi_ror2_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$ror2$obs);

## Make a helper function to set cline models
tg.loadror2model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$ror2$models[[id]] <<- hzar.makeCline1DFreq(tg$ror2$obs, scaling, tails)

tg.loadror2model("fixed","none","modelI");
tg.loadror2model("free" ,"none","modelII");
tg.loadror2model("free" ,"both","modelIII");
tg.loadror2model("free" ,"right","modelIV");
tg.loadror2model("free" ,"left","modelV");
tg.loadror2model("free" ,"mirror","modelVI");

## Check the default settings
print(tg$ror2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$ror2$models <- sapply(tg$ror2$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$ror2$models)

## Compile each of the models to prepare for fitting
tg$ror2$fitRs$init <- sapply(tg$ror2$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$ror2$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$ror2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ror2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ror2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$ror2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ror2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ror2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$ror2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ror2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ror2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$ror2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ror2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ror2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$ror2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ror2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ror2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$ror2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ror2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ror2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$ror2$fitRs$init)

## Run the models for an initial chain
tg$ror2$runs$init <- list()

tg$ror2$runs$init$modelI <-
  hzar.doFit(tg$ror2$fitRs$init$modelI)

tg$ror2$runs$init$modelII <-
  hzar.doFit(tg$ror2$fitRs$init$modelII)

tg$ror2$runs$init$modelIII <-
  hzar.doFit(tg$ror2$fitRs$init$modelIII)

tg$ror2$runs$init$modelIV <-
  hzar.doFit(tg$ror2$fitRs$init$modelIV)

tg$ror2$runs$init$modelV <-
  hzar.doFit(tg$ror2$fitRs$init$modelV)

tg$ror2$runs$init$modelVI <-
  hzar.doFit(tg$ror2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$ror2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$ror2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$ror2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$ror2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$ror2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$ror2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$ror2$fitRs$chains <-
  lapply(tg$ror2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$ror2$fitRs$chains <-
  hzar.multiFitRequest(tg$ror2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$ror2$runs$chains <-  hzar.doChain.multi(tg$ror2$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$ror2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$ror2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$ror2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$ror2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$ror2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$ror2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$ror2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$ror2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$ror2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$ror2$runs$init$modelI)
tg$ror2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$ror2$runs$init$modelII)
tg$ror2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$ror2$runs$init$modelIII)
tg$ror2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$ror2$runs$init$modelIV)
tg$ror2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$ror2$runs$init$modelV)
tg$ror2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$ror2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$ror2$analysis$oDG <-
  hzar.make.obsDataGroup(tg$ror2$analysis$initDGs)
tg$ror2$analysis$oDG <-
  hzar.copyModelLabels(tg$ror2$analysis$initDGs,
                       tg$ror2$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$ror2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$ror2$runs$chains,
                                hzar.dataGroup.add),
                         tg$ror2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$ror2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$ror2$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$ror2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$ror2$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$ror2$analysis$model.name <-
        rownames(tg$ror2$analysis$AICcTable
        )[[ which.min(tg$ror2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$ror2$analysis$model.selected <-
  tg$ror2$analysis$oDG$data.groups[[tg$ror2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$ror2$analysis$model.selected,
                         names(tg$ror2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$ror2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$ror2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$ror2$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on APC-CAMK4 region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$apc.camk4 <- list();
## Space to hold the observed data
tg$apc.camk4$obs <- list();
## Space to hold the models to fit
tg$apc.camk4$models <- list();
## Space to hold the compiled fit requests
tg$apc.camk4$fitRs <- list();
## Space to hold the output data chains
tg$apc.camk4$runs <- list();
## Space to hold the analysed data
tg$apc.camk4$analysis <- list();

## Assign data for apc.camk4
tg$apc.camk4$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_apc_camk4_mean,
                             tg.dists$hi_apc_camk4_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$apc.camk4$obs);

## Make a helper function to set cline models
tg.loadapc.camk4model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$apc.camk4$models[[id]] <<- hzar.makeCline1DFreq(tg$apc.camk4$obs, scaling, tails)

tg.loadapc.camk4model("fixed","none","modelI");
tg.loadapc.camk4model("free" ,"none","modelII");
tg.loadapc.camk4model("free" ,"both","modelIII");
tg.loadapc.camk4model("free" ,"right","modelIV");
tg.loadapc.camk4model("free" ,"left","modelV");
tg.loadapc.camk4model("free" ,"mirror","modelVI");

## Check the default settings
print(tg$apc.camk4$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$apc.camk4$models <- sapply(tg$apc.camk4$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$apc.camk4$models)

## Compile each of the models to prepare for fitting
tg$apc.camk4$fitRs$init <- sapply(tg$apc.camk4$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$apc.camk4$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$apc.camk4$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$apc.camk4$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$apc.camk4$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$apc.camk4$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$apc.camk4$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$apc.camk4$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$apc.camk4$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$apc.camk4$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$apc.camk4$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$apc.camk4$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$apc.camk4$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$apc.camk4$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$apc.camk4$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$apc.camk4$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$apc.camk4$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$apc.camk4$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$apc.camk4$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$apc.camk4$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$apc.camk4$fitRs$init)

## Run the models for an initial chain
tg$apc.camk4$runs$init <- list()

tg$apc.camk4$runs$init$modelI <-
  hzar.doFit(tg$apc.camk4$fitRs$init$modelI)

tg$apc.camk4$runs$init$modelII <-
  hzar.doFit(tg$apc.camk4$fitRs$init$modelII)

tg$apc.camk4$runs$init$modelIII <-
  hzar.doFit(tg$apc.camk4$fitRs$init$modelIII)

tg$apc.camk4$runs$init$modelIV <-
  hzar.doFit(tg$apc.camk4$fitRs$init$modelIV)

tg$apc.camk4$runs$init$modelV <-
  hzar.doFit(tg$apc.camk4$fitRs$init$modelV)

tg$apc.camk4$runs$init$modelVI <-
  hzar.doFit(tg$apc.camk4$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$apc.camk4$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$apc.camk4$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$apc.camk4$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$apc.camk4$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$apc.camk4$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$apc.camk4$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$apc.camk4$fitRs$chains <-
  lapply(tg$apc.camk4$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$apc.camk4$fitRs$chains <-
  hzar.multiFitRequest(tg$apc.camk4$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$apc.camk4$runs$chains <-  hzar.doChain.multi(tg$apc.camk4$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$apc.camk4$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$apc.camk4$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$apc.camk4$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$apc.camk4$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$apc.camk4$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$apc.camk4$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$apc.camk4$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$apc.camk4$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$apc.camk4$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$apc.camk4$runs$init$modelI)
tg$apc.camk4$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$apc.camk4$runs$init$modelII)
tg$apc.camk4$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$apc.camk4$runs$init$modelIII)
tg$apc.camk4$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$apc.camk4$runs$init$modelIV)
tg$apc.camk4$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$apc.camk4$runs$init$modelV)
tg$apc.camk4$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$apc.camk4$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$apc.camk4$analysis$oDG <-
  hzar.make.obsDataGroup(tg$apc.camk4$analysis$initDGs)
tg$apc.camk4$analysis$oDG <-
  hzar.copyModelLabels(tg$apc.camk4$analysis$initDGs,
                       tg$apc.camk4$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$apc.camk4$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$apc.camk4$runs$chains,
                                hzar.dataGroup.add),
                         tg$apc.camk4$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$apc.camk4$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$apc.camk4$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$apc.camk4$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$apc.camk4$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$apc.camk4$analysis$model.name <-
        rownames(tg$apc.camk4$analysis$AICcTable
        )[[ which.min(tg$apc.camk4$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$apc.camk4$analysis$model.selected <-
  tg$apc.camk4$analysis$oDG$data.groups[[tg$apc.camk4$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$apc.camk4$analysis$model.selected,
                         names(tg$apc.camk4$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$apc.camk4$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$apc.camk4$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$apc.camk4$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on ICE1-lncRNA region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$ice1.lncrna <- list();
## Space to hold the observed data
tg$ice1.lncrna$obs <- list();
## Space to hold the models to fit
tg$ice1.lncrna$models <- list();
## Space to hold the compiled fit requests
tg$ice1.lncrna$fitRs <- list();
## Space to hold the output data chains
tg$ice1.lncrna$runs <- list();
## Space to hold the analysed data
tg$ice1.lncrna$analysis <- list();

## Assign data for ice1.lncrna
tg$ice1.lncrna$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_ice1_lncrna_mean,
                             tg.dists$hi_ice1_lncrna_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$ice1.lncrna$obs);

## Make a helper function to set cline models
tg.loadice1.lncrnamodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$ice1.lncrna$models[[id]] <<- hzar.makeCline1DFreq(tg$ice1.lncrna$obs, scaling, tails)

tg.loadice1.lncrnamodel("fixed","none","modelI");
tg.loadice1.lncrnamodel("free" ,"none","modelII");
tg.loadice1.lncrnamodel("free" ,"both","modelIII");
tg.loadice1.lncrnamodel("free" ,"right","modelIV");
tg.loadice1.lncrnamodel("free" ,"left","modelV");
tg.loadice1.lncrnamodel("free" ,"mirror","modelVI");

## Check the default settings
print(tg$ice1.lncrna$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$ice1.lncrna$models <- sapply(tg$ice1.lncrna$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$ice1.lncrna$models)

## Compile each of the models to prepare for fitting
tg$ice1.lncrna$fitRs$init <- sapply(tg$ice1.lncrna$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$ice1.lncrna$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$ice1.lncrna$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ice1.lncrna$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ice1.lncrna$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$ice1.lncrna$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ice1.lncrna$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ice1.lncrna$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$ice1.lncrna$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ice1.lncrna$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ice1.lncrna$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$ice1.lncrna$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ice1.lncrna$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ice1.lncrna$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$ice1.lncrna$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ice1.lncrna$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ice1.lncrna$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$ice1.lncrna$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$ice1.lncrna$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$ice1.lncrna$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$ice1.lncrna$fitRs$init)

## Run the models for an initial chain
tg$ice1.lncrna$runs$init <- list()

tg$ice1.lncrna$runs$init$modelI <-
  hzar.doFit(tg$ice1.lncrna$fitRs$init$modelI)

tg$ice1.lncrna$runs$init$modelII <-
  hzar.doFit(tg$ice1.lncrna$fitRs$init$modelII)

tg$ice1.lncrna$runs$init$modelIII <-
  hzar.doFit(tg$ice1.lncrna$fitRs$init$modelIII)

tg$ice1.lncrna$runs$init$modelIV <-
  hzar.doFit(tg$ice1.lncrna$fitRs$init$modelIV)

tg$ice1.lncrna$runs$init$modelV <-
  hzar.doFit(tg$ice1.lncrna$fitRs$init$modelV)

tg$ice1.lncrna$runs$init$modelVI <-
  hzar.doFit(tg$ice1.lncrna$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$ice1.lncrna$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$ice1.lncrna$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$ice1.lncrna$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$ice1.lncrna$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$ice1.lncrna$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$ice1.lncrna$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$ice1.lncrna$fitRs$chains <-
  lapply(tg$ice1.lncrna$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$ice1.lncrna$fitRs$chains <-
  hzar.multiFitRequest(tg$ice1.lncrna$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$ice1.lncrna$runs$chains <-  hzar.doChain.multi(tg$ice1.lncrna$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$ice1.lncrna$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$ice1.lncrna$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$ice1.lncrna$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$ice1.lncrna$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$ice1.lncrna$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$ice1.lncrna$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$ice1.lncrna$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$ice1.lncrna$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$ice1.lncrna$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$ice1.lncrna$runs$init$modelI)
tg$ice1.lncrna$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$ice1.lncrna$runs$init$modelII)
tg$ice1.lncrna$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$ice1.lncrna$runs$init$modelIII)
tg$ice1.lncrna$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$ice1.lncrna$runs$init$modelIV)
tg$ice1.lncrna$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$ice1.lncrna$runs$init$modelV)
tg$ice1.lncrna$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$ice1.lncrna$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$ice1.lncrna$analysis$oDG <-
  hzar.make.obsDataGroup(tg$ice1.lncrna$analysis$initDGs)
tg$ice1.lncrna$analysis$oDG <-
  hzar.copyModelLabels(tg$ice1.lncrna$analysis$initDGs,
                       tg$ice1.lncrna$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$ice1.lncrna$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$ice1.lncrna$runs$chains,
                                hzar.dataGroup.add),
                         tg$ice1.lncrna$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$ice1.lncrna$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$ice1.lncrna$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$ice1.lncrna$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$ice1.lncrna$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$ice1.lncrna$analysis$model.name <-
        rownames(tg$ice1.lncrna$analysis$AICcTable
        )[[ which.min(tg$ice1.lncrna$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$ice1.lncrna$analysis$model.selected <-
  tg$ice1.lncrna$analysis$oDG$data.groups[[tg$ice1.lncrna$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$ice1.lncrna$analysis$model.selected,
                         names(tg$ice1.lncrna$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$ice1.lncrna$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$ice1.lncrna$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$ice1.lncrna$analysis$model.selected);

## End Analysis

#dev.off()


### Analysis on PDE1C region (rustica-tytleri)------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg$pde1c <- list();
## Space to hold the observed data
tg$pde1c$obs <- list();
## Space to hold the models to fit
tg$pde1c$models <- list();
## Space to hold the compiled fit requests
tg$pde1c$fitRs <- list();
## Space to hold the output data chains
tg$pde1c$runs <- list();
## Space to hold the analysed data
tg$pde1c$analysis <- list();

## Assign data for pde1c
tg$pde1c$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$hi_pde1c_mean,
                             tg.dists$hi_pde1c_count);

## Look at a graph of the observed data
hzar.plot.obsData(tg$pde1c$obs);

## Make a helper function to set cline models
tg.loadpde1cmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  tg$pde1c$models[[id]] <<- hzar.makeCline1DFreq(tg$pde1c$obs, scaling, tails)

tg.loadpde1cmodel("fixed","none","modelI");
tg.loadpde1cmodel("free" ,"none","modelII");
tg.loadpde1cmodel("free" ,"both","modelIII");
tg.loadpde1cmodel("free" ,"right","modelIV");
tg.loadpde1cmodel("free" ,"left","modelV");
tg.loadpde1cmodel("free" ,"mirror","modelVI");

## Check the default settings
print(tg$pde1c$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(tg.dists$km)
max(tg.dists$km)

tg$pde1c$models <- sapply(tg$pde1c$models,
                          hzar.model.addBoxReq,
                          -30 , 2185,
                          simplify=FALSE)

## Check the updated settings
print(tg$pde1c$models)

## Compile each of the models to prepare for fitting
tg$pde1c$fitRs$init <- sapply(tg$pde1c$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=tg$pde1c$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
tg$pde1c$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$pde1c$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$pde1c$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg$pde1c$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$pde1c$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$pde1c$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg$pde1c$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$pde1c$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$pde1c$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$pde1c$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$pde1c$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$pde1c$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$pde1c$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$pde1c$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$pde1c$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg$pde1c$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg$pde1c$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg$pde1c$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg$pde1c$fitRs$init)

## Run the models for an initial chain
tg$pde1c$runs$init <- list()

tg$pde1c$runs$init$modelI <-
  hzar.doFit(tg$pde1c$fitRs$init$modelI)

tg$pde1c$runs$init$modelII <-
  hzar.doFit(tg$pde1c$fitRs$init$modelII)

tg$pde1c$runs$init$modelIII <-
  hzar.doFit(tg$pde1c$fitRs$init$modelIII)

tg$pde1c$runs$init$modelIV <-
  hzar.doFit(tg$pde1c$fitRs$init$modelIV)

tg$pde1c$runs$init$modelV <-
  hzar.doFit(tg$pde1c$fitRs$init$modelV)

tg$pde1c$runs$init$modelVI <-
  hzar.doFit(tg$pde1c$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(tg$pde1c$runs$init$modelI))
#plot(hzar.mcmc.bindLL(tg$pde1c$runs$init$modelII))
#plot(hzar.mcmc.bindLL(tg$pde1c$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(tg$pde1c$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(tg$pde1c$runs$init$modelV))
#plot(hzar.mcmc.bindLL(tg$pde1c$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
tg$pde1c$fitRs$chains <-
  lapply(tg$pde1c$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg$pde1c$fitRs$chains <-
  hzar.multiFitRequest(tg$pde1c$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg$pde1c$runs$chains <-  hzar.doChain.multi(tg$pde1c$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg$pde1c$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg$pde1c$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg$pde1c$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg$pde1c$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg$pde1c$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg$pde1c$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg$pde1c$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg$pde1c$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg$pde1c$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg$pde1c$runs$init$modelI)
tg$pde1c$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg$pde1c$runs$init$modelII)
tg$pde1c$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg$pde1c$runs$init$modelIII)
tg$pde1c$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg$pde1c$runs$init$modelIV)
tg$pde1c$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg$pde1c$runs$init$modelV)
tg$pde1c$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg$pde1c$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg$pde1c$analysis$oDG <-
  hzar.make.obsDataGroup(tg$pde1c$analysis$initDGs)
tg$pde1c$analysis$oDG <-
  hzar.copyModelLabels(tg$pde1c$analysis$initDGs,
                       tg$pde1c$analysis$oDG)

## Convetg all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg$pde1c$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg$pde1c$runs$chains,
                                hzar.dataGroup.add),
                         tg$pde1c$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg$pde1c$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(tg$pde1c$analysis$oDG);

## Do model selection based on the AICc scores
print(tg$pde1c$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$pde1c$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg$pde1c$analysis$model.name <-
        rownames(tg$pde1c$analysis$AICcTable
        )[[ which.min(tg$pde1c$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg$pde1c$analysis$model.selected <-
  tg$pde1c$analysis$oDG$data.groups[[tg$pde1c$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$pde1c$analysis$model.selected,
                         names(tg$pde1c$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$pde1c$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$pde1c$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(tg$pde1c$analysis$model.selected);

## End Analysis

#dev.off()

### Save analysis results to Rdata object-----------------------------------

save(tg, file = './candidate_results_hzar/hzar_data_tytleri-gutturalis.RData')
