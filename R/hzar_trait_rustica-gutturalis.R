############################################################################
# Barn swallow geographic clines in hybrid zone transects
############################################################################

# This script contains commands for analysis of individual hybrid index and
# geographic clines for 'outlier' regions of the genome versus background
# to understand if genomic regions associated with mate choice traits show
# restricted gene flow (i.e., steeper/narrower clines).

# This workflow involves:
# 1. Estimation of hybrid indices from outlier regions for the rustica-gutturalis hybrid zone
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
load("./candidate_results_hzar/hzar_data_rustica-gutturalis.RData")

### Introgress - estimation of hybrid index (rustica-gutturalis)------------

## Read in VCF data

rg_ru.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.rustica-gutturalis_rustica.vcf.gz')
rg_gu.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.rustica-gutturalis_rustica.vcf.gz')
rg_gu.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.rustica-gutturalis_rustica.vcf.gz')
rg_gu.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.rustica-gutturalis_rustica.vcf.gz')
rg_gu.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.rustica-gutturalis_rustica.vcf.gz')
rg_gu.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.rustica-gutturalis_rustica.vcf.gz')
rg_gu.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.rustica-gutturalis_rustica.vcf.gz')
rg_gu.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.rustica-gutturalis_rustica.vcf.gz')
rg_gu.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.rustica-gutturalis_rustica.vcf.gz')
rg_gu.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.rustica-gutturalis_hybrids.vcf.gz')

rg_ru.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.rustica-gutturalis_rustica.vcf.gz')
rg_gu.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.rustica-gutturalis_gutturalis.vcf.gz')
rg_rg.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.rustica-gutturalis_hybrids.vcf.gz')

## Extract genotype matrices

# Outlier regions
rg_ru.kitlg.gt <- extract.gt(rg_ru.kitlg.vcfr)
rg_gu.kitlg.gt <- extract.gt(rg_gu.kitlg.vcfr)
rg_rg.kitlg.gt <- extract.gt(rg_rg.kitlg.vcfr)

rg_ru.plxnc1.gt <- extract.gt(rg_ru.plxnc1.vcfr)
rg_gu.plxnc1.gt <- extract.gt(rg_gu.plxnc1.vcfr)
rg_rg.plxnc1.gt <- extract.gt(rg_rg.plxnc1.vcfr)

rg_ru.spef2.prlr.gt <- extract.gt(rg_ru.spef2.prlr.vcfr)
rg_gu.spef2.prlr.gt <- extract.gt(rg_gu.spef2.prlr.vcfr)
rg_rg.spef2.prlr.gt <- extract.gt(rg_rg.spef2.prlr.vcfr)

rg_ru.slc45a2.gt <- extract.gt(rg_ru.slc45a2.vcfr)
rg_gu.slc45a2.gt <- extract.gt(rg_gu.slc45a2.vcfr)
rg_rg.slc45a2.gt <- extract.gt(rg_rg.slc45a2.vcfr)

rg_ru.bnc2.gt <- extract.gt(rg_ru.bnc2.vcfr)
rg_gu.bnc2.gt <- extract.gt(rg_gu.bnc2.vcfr)
rg_rg.bnc2.gt <- extract.gt(rg_rg.bnc2.vcfr)

rg_ru.gnaq.gt <- extract.gt(rg_ru.gnaq.vcfr)
rg_gu.gnaq.gt <- extract.gt(rg_gu.gnaq.vcfr)
rg_rg.gnaq.gt <- extract.gt(rg_rg.gnaq.vcfr)

rg_ru.ror2.gt <- extract.gt(rg_ru.ror2.vcfr)
rg_gu.ror2.gt <- extract.gt(rg_gu.ror2.vcfr)
rg_rg.ror2.gt <- extract.gt(rg_rg.ror2.vcfr)

rg_ru.apc.camk4.gt <- extract.gt(rg_ru.apc.camk4.vcfr)
rg_gu.apc.camk4.gt <- extract.gt(rg_gu.apc.camk4.vcfr)
rg_rg.apc.camk4.gt <- extract.gt(rg_rg.apc.camk4.vcfr)

rg_ru.ice1.lncrna.gt <- extract.gt(rg_ru.ice1.lncrna.vcfr)
rg_gu.ice1.lncrna.gt <- extract.gt(rg_gu.ice1.lncrna.vcfr)
rg_rg.ice1.lncrna.gt <- extract.gt(rg_rg.ice1.lncrna.vcfr)

rg_ru.pde1c.gt <- extract.gt(rg_ru.pde1c.vcfr)
rg_gu.pde1c.gt <- extract.gt(rg_gu.pde1c.vcfr)
rg_rg.pde1c.gt <- extract.gt(rg_rg.pde1c.vcfr)

## Set locus sets

rg.kitlg.loci.dat <- as.data.frame(rg_rg.kitlg.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.plxnc1.loci.dat <- as.data.frame(rg_rg.plxnc1.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.spef2.prlr.loci.dat <- as.data.frame(rg_rg.spef2.prlr.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.slc45a2.loci.dat <- as.data.frame(rg_rg.slc45a2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.bnc2.loci.dat <- as.data.frame(rg_rg.bnc2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.gnaq.loci.dat <- as.data.frame(rg_rg.gnaq.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.ror2.loci.dat <- as.data.frame(rg_rg.ror2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.apc.camk4.loci.dat <- as.data.frame(rg_rg.apc.camk4.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.ice1.lncrna.loci.dat <- as.data.frame(rg_rg.ice1.lncrna.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rg.pde1c.loci.dat <- as.data.frame(rg_rg.pde1c.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

## Prepare data for analysis

rg_ru.kitlg.dat <- prepare.data(rg_ru.kitlg.gt, loci.data = rg.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.kitlg.gt, parental2 = rg_gu.kitlg.gt)
rg_gu.kitlg.dat <- prepare.data(rg_gu.kitlg.gt, loci.data = rg.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.kitlg.gt, parental2 = rg_gu.kitlg.gt)
rg_rg.kitlg.dat <- prepare.data(rg_rg.kitlg.gt, loci.data = rg.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.kitlg.gt, parental2 = rg_gu.kitlg.gt)

rg_ru.plxnc1.dat <- prepare.data(rg_ru.plxnc1.gt, loci.data = rg.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.plxnc1.gt, parental2 = rg_gu.plxnc1.gt)
rg_gu.plxnc1.dat <- prepare.data(rg_gu.plxnc1.gt, loci.data = rg.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.plxnc1.gt, parental2 = rg_gu.plxnc1.gt)
rg_rg.plxnc1.dat <- prepare.data(rg_rg.plxnc1.gt, loci.data = rg.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.plxnc1.gt, parental2 = rg_gu.plxnc1.gt)

rg_ru.spef2.prlr.dat <- prepare.data(rg_ru.spef2.prlr.gt, loci.data = rg.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.spef2.prlr.gt, parental2 = rg_gu.spef2.prlr.gt)
rg_gu.spef2.prlr.dat <- prepare.data(rg_gu.spef2.prlr.gt, loci.data = rg.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.spef2.prlr.gt, parental2 = rg_gu.spef2.prlr.gt)
rg_rg.spef2.prlr.dat <- prepare.data(rg_rg.spef2.prlr.gt, loci.data = rg.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.spef2.prlr.gt, parental2 = rg_gu.spef2.prlr.gt)

rg_ru.slc45a2.dat <- prepare.data(rg_ru.slc45a2.gt, loci.data = rg.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.slc45a2.gt, parental2 = rg_gu.slc45a2.gt)
rg_gu.slc45a2.dat <- prepare.data(rg_gu.slc45a2.gt, loci.data = rg.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.slc45a2.gt, parental2 = rg_gu.slc45a2.gt)
rg_rg.slc45a2.dat <- prepare.data(rg_rg.slc45a2.gt, loci.data = rg.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.slc45a2.gt, parental2 = rg_gu.slc45a2.gt)

rg_ru.bnc2.dat <- prepare.data(rg_ru.bnc2.gt, loci.data = rg.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.bnc2.gt, parental2 = rg_gu.bnc2.gt)
rg_gu.bnc2.dat <- prepare.data(rg_gu.bnc2.gt, loci.data = rg.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.bnc2.gt, parental2 = rg_gu.bnc2.gt)
rg_rg.bnc2.dat <- prepare.data(rg_rg.bnc2.gt, loci.data = rg.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.bnc2.gt, parental2 = rg_gu.bnc2.gt)

rg_ru.gnaq.dat <- prepare.data(rg_ru.gnaq.gt, loci.data = rg.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.gnaq.gt, parental2 = rg_gu.gnaq.gt)
rg_gu.gnaq.dat <- prepare.data(rg_gu.gnaq.gt, loci.data = rg.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.gnaq.gt, parental2 = rg_gu.gnaq.gt)
rg_rg.gnaq.dat <- prepare.data(rg_rg.gnaq.gt, loci.data = rg.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.gnaq.gt, parental2 = rg_gu.gnaq.gt)

rg_ru.ror2.dat <- prepare.data(rg_ru.ror2.gt, loci.data = rg.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.ror2.gt, parental2 = rg_gu.ror2.gt)
rg_gu.ror2.dat <- prepare.data(rg_gu.ror2.gt, loci.data = rg.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.ror2.gt, parental2 = rg_gu.ror2.gt)
rg_rg.ror2.dat <- prepare.data(rg_rg.ror2.gt, loci.data = rg.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.ror2.gt, parental2 = rg_gu.ror2.gt)

rg_ru.apc.camk4.dat <- prepare.data(rg_ru.apc.camk4.gt, loci.data = rg.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.apc.camk4.gt, parental2 = rg_gu.apc.camk4.gt)
rg_gu.apc.camk4.dat <- prepare.data(rg_gu.apc.camk4.gt, loci.data = rg.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.apc.camk4.gt, parental2 = rg_gu.apc.camk4.gt)
rg_rg.apc.camk4.dat <- prepare.data(rg_rg.apc.camk4.gt, loci.data = rg.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.apc.camk4.gt, parental2 = rg_gu.apc.camk4.gt)

rg_ru.ice1.lncrna.dat <- prepare.data(rg_ru.ice1.lncrna.gt, loci.data = rg.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.ice1.lncrna.gt, parental2 = rg_gu.ice1.lncrna.gt)
rg_gu.ice1.lncrna.dat <- prepare.data(rg_gu.ice1.lncrna.gt, loci.data = rg.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.ice1.lncrna.gt, parental2 = rg_gu.ice1.lncrna.gt)
rg_rg.ice1.lncrna.dat <- prepare.data(rg_rg.ice1.lncrna.gt, loci.data = rg.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.ice1.lncrna.gt, parental2 = rg_gu.ice1.lncrna.gt)

rg_ru.pde1c.dat <- prepare.data(rg_ru.pde1c.gt, loci.data = rg.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.pde1c.gt, parental2 = rg_gu.pde1c.gt)
rg_gu.pde1c.dat <- prepare.data(rg_gu.pde1c.gt, loci.data = rg.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.pde1c.gt, parental2 = rg_gu.pde1c.gt)
rg_rg.pde1c.dat <- prepare.data(rg_rg.pde1c.gt, loci.data = rg.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rg_ru.pde1c.gt, parental2 = rg_gu.pde1c.gt)

## Estimate hybrid index

rg_ru.kitlg.h <- est.h(rg_ru.kitlg.dat, loci.data = rg.kitlg.loci.dat, fixed = FALSE)
rg_gu.kitlg.h <- est.h(rg_gu.kitlg.dat, loci.data = rg.kitlg.loci.dat, fixed = FALSE)
rg_rg.kitlg.h <- est.h(rg_rg.kitlg.dat, loci.data = rg.kitlg.loci.dat, fixed = FALSE)

rg_ru.plxnc1.h <- est.h(rg_ru.plxnc1.dat, loci.data = rg.plxnc1.loci.dat, fixed = FALSE)
rg_gu.plxnc1.h <- est.h(rg_gu.plxnc1.dat, loci.data = rg.plxnc1.loci.dat, fixed = FALSE)
rg_rg.plxnc1.h <- est.h(rg_rg.plxnc1.dat, loci.data = rg.plxnc1.loci.dat, fixed = FALSE)

rg_ru.spef2.prlr.h <- est.h(rg_ru.spef2.prlr.dat, loci.data = rg.spef2.prlr.loci.dat, fixed = FALSE)
rg_gu.spef2.prlr.h <- est.h(rg_gu.spef2.prlr.dat, loci.data = rg.spef2.prlr.loci.dat, fixed = FALSE)
rg_rg.spef2.prlr.h <- est.h(rg_rg.spef2.prlr.dat, loci.data = rg.spef2.prlr.loci.dat, fixed = FALSE)

rg_ru.slc45a2.h <- est.h(rg_ru.slc45a2.dat, loci.data = rg.slc45a2.loci.dat, fixed = FALSE)
rg_gu.slc45a2.h <- est.h(rg_gu.slc45a2.dat, loci.data = rg.slc45a2.loci.dat, fixed = FALSE)
rg_rg.slc45a2.h <- est.h(rg_rg.slc45a2.dat, loci.data = rg.slc45a2.loci.dat, fixed = FALSE)

rg_ru.bnc2.h <- est.h(rg_ru.bnc2.dat, loci.data = rg.bnc2.loci.dat, fixed = FALSE)
rg_gu.bnc2.h <- est.h(rg_gu.bnc2.dat, loci.data = rg.bnc2.loci.dat, fixed = FALSE)
rg_rg.bnc2.h <- est.h(rg_rg.bnc2.dat, loci.data = rg.bnc2.loci.dat, fixed = FALSE)

rg_ru.gnaq.h <- est.h(rg_ru.gnaq.dat, loci.data = rg.gnaq.loci.dat, fixed = FALSE)
rg_gu.gnaq.h <- est.h(rg_gu.gnaq.dat, loci.data = rg.gnaq.loci.dat, fixed = FALSE)
rg_rg.gnaq.h <- est.h(rg_rg.gnaq.dat, loci.data = rg.gnaq.loci.dat, fixed = FALSE)

rg_ru.ror2.h <- est.h(rg_ru.ror2.dat, loci.data = rg.ror2.loci.dat, fixed = FALSE)
rg_gu.ror2.h <- est.h(rg_gu.ror2.dat, loci.data = rg.ror2.loci.dat, fixed = FALSE)
rg_rg.ror2.h <- est.h(rg_rg.ror2.dat, loci.data = rg.ror2.loci.dat, fixed = FALSE)

rg_ru.apc.camk4.h <- est.h(rg_ru.apc.camk4.dat, loci.data = rg.apc.camk4.loci.dat, fixed = FALSE)
rg_gu.apc.camk4.h <- est.h(rg_gu.apc.camk4.dat, loci.data = rg.apc.camk4.loci.dat, fixed = FALSE)
rg_rg.apc.camk4.h <- est.h(rg_rg.apc.camk4.dat, loci.data = rg.apc.camk4.loci.dat, fixed = FALSE)

rg_ru.ice1.lncrna.h <- est.h(rg_ru.ice1.lncrna.dat, loci.data = rg.ice1.lncrna.loci.dat, fixed = FALSE)
rg_gu.ice1.lncrna.h <- est.h(rg_gu.ice1.lncrna.dat, loci.data = rg.ice1.lncrna.loci.dat, fixed = FALSE)
rg_rg.ice1.lncrna.h <- est.h(rg_rg.ice1.lncrna.dat, loci.data = rg.ice1.lncrna.loci.dat, fixed = FALSE)

rg_ru.pde1c.h <- est.h(rg_ru.pde1c.dat, loci.data = rg.pde1c.loci.dat, fixed = FALSE)
rg_gu.pde1c.h <- est.h(rg_gu.pde1c.dat, loci.data = rg.pde1c.loci.dat, fixed = FALSE)
rg_rg.pde1c.h <- est.h(rg_rg.pde1c.dat, loci.data = rg.pde1c.loci.dat, fixed = FALSE)

## Query hybrid index values and write to `data.csv`

print(rg_ru.kitlg.h)
print(rg_rg.kitlg.h)
print(rg_gu.kitlg.h)

print(rg_ru.plxnc1.h)
print(rg_rg.plxnc1.h)
print(rg_gu.plxnc1.h)

print(rg_ru.spef2.prlr.h)
print(rg_rg.spef2.prlr.h)
print(rg_gu.spef2.prlr.h)

print(rg_ru.slc45a2.h)
print(rg_rg.slc45a2.h)
print(rg_gu.slc45a2.h)

print(rg_ru.bnc2.h)
print(rg_rg.bnc2.h)
print(rg_gu.bnc2.h)

print(rg_ru.gnaq.h)
print(rg_rg.gnaq.h)
print(rg_gu.gnaq.h)

print(rg_ru.ror2.h)
print(rg_rg.ror2.h)
print(rg_gu.ror2.h)

print(rg_ru.apc.camk4.h)
print(rg_rg.apc.camk4.h)
print(rg_gu.apc.camk4.h)

print(rg_ru.ice1.lncrna.h)
print(rg_rg.ice1.lncrna.h)
print(rg_gu.ice1.lncrna.h)

print(rg_ru.pde1c.h)
print(rg_rg.pde1c.h)
print(rg_gu.pde1c.h)

### Read in and format data for HZAR (rustica-gutturalis)--------------------

# These data include the sample, population, lat, long, locality, hybrid index
# (background), and hybrid index (outlier regions).

data.rg <- read.csv('./data.rustica-gutturalis.csv',header=T)
names(data.rg)

## Standardize phenotype vector ranges between 0 and 1 and add to data frame
mean_rwl_norm = (data.rg$mean_rwl-min(data.rg$mean_rwl,na.rm = T))/(max(data.rg$mean_rwl,na.rm = T)-min(data.rg$mean_rwl,na.rm = T))
mean_ts_random_norm = (data.rg$mean_ts_random-min(data.rg$mean_ts_random,na.rm = T))/(max(data.rg$mean_ts_random,na.rm = T)-min(data.rg$mean_ts_random,na.rm = T))
breast_avg_bright_norm = (data.rg$breast_avg_bright-min(data.rg$breast_avg_bright,na.rm = T))/(max(data.rg$breast_avg_bright,na.rm = T)-min(data.rg$breast_avg_bright,na.rm = T))

data.rg$wing <- mean_rwl_norm
data.rg$tail <- mean_ts_random_norm
data.rg$vent <- breast_avg_bright_norm

## Subset to get locality, long, lat, HI-background, and HI-outlier
data.rg <- dplyr::select(data.rg, loc,long, lat, hi_kitlg, hi_plxnc1, hi_spef2_prlr, hi_slc45a2, hi_bnc2, hi_gnaq, hi_ror2, hi_apc_camk4, hi_ice1_lncrna, hi_pde1c)
head(data.rg)

## Calculate means, variances, and counts

data.rg.melt<-reshape2::melt(data.rg, id.vars=c("loc"))
str(data.rg.melt)
head(data.rg.melt)
data.rg.hzar<- data.rg.melt %>%
  group_by(loc, variable) %>%
  dplyr::summarise(mean=mean(value, na.rm=TRUE), var=var(value, na.rm=TRUE), count=sum(!is.na(value)))
head(data.rg.hzar)

## Change to wide format

d<-reshape2::melt(data.rg.hzar, id.vars=c("loc", "variable"))
head(d)
colnames(d)[3]<-"metric"
data.rg.hzar.cast<-dcast(d, loc~ variable + metric )

print(data.rg.hzar.cast)

## Calculate distances (km) for transect

rg.dists<-hzar.map.latLongSites(data.rg.hzar.cast$loc, data.rg.hzar.cast$lat_mean, data.rg.hzar.cast$long_mean, degrees = TRUE)

join_func<-function(y){
  y<-left_join(y, data.rg.hzar.cast, by=c("site"="loc"))
  y<-y[,-c(2:3)]
  y<-arrange(y, km)
}

## Subset points to the sampling transects

rg.dists<-filter(rg.dists, site=="karasuk" | site=="urumqi" | site=="wuwei" | site=="Lanzhou" | site=="Jiuquan" | site=="Zhangye" | site=="gaotai" | site=="yumen" | site=="xian" | site=="zhengzhou") #urumqi the straight line
print(rg.dists)

##Add km transect distances

rg.dists$km <- hzar.map.distanceFromSite(rg.dists,"karasuk",units="Km")
rg.dists<-join_func(rg.dists)
#rg.dists<-na.omit(rg.dists) # To keep urumqi in there. UNCOMMENT IF THIS MESSES UP ANCESTRY-BASED ANALYSES!!!
print(rg.dists)

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

### Analysis on KITLG region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$kitlg <- list();
## Space to hold the observed data
rg$kitlg$obs <- list();
## Space to hold the models to fit
rg$kitlg$models <- list();
## Space to hold the compiled fit requests
rg$kitlg$fitRs <- list();
## Space to hold the output data chains
rg$kitlg$runs <- list();
## Space to hold the analysed data
rg$kitlg$analysis <- list();

## Assign data for KITLG region
rg$kitlg$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_kitlg_mean,
                             rg.dists$hi_kitlg_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$kitlg$obs);

## Make a helper function to set cline models
rg.loadkitlgmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$kitlg$models[[id]] <<- hzar.makeCline1DFreq(rg$kitlg$obs, scaling, tails)

rg.loadkitlgmodel("fixed","none","modelI");
rg.loadkitlgmodel("free" ,"none","modelII");
rg.loadkitlgmodel("free" ,"both","modelIII");
rg.loadkitlgmodel("free" ,"right","modelIV");
rg.loadkitlgmodel("free" ,"left","modelV");
rg.loadkitlgmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rg$kitlg$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$kitlg$models <- sapply(rg$kitlg$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$kitlg$models)

## Compile each of the models to prepare for fitting
rg$kitlg$fitRs$init <- sapply(rg$kitlg$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$kitlg$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$kitlg$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$kitlg$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$kitlg$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$kitlg$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$kitlg$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$kitlg$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$kitlg$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$kitlg$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$kitlg$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$kitlg$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$kitlg$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$kitlg$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$kitlg$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$kitlg$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$kitlg$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$kitlg$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$kitlg$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$kitlg$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$kitlg$fitRs$init)

## Run the models for an initial chain
rg$kitlg$runs$init <- list()

rg$kitlg$runs$init$modelI <-
  hzar.doFit(rg$kitlg$fitRs$init$modelI)

rg$kitlg$runs$init$modelII <-
  hzar.doFit(rg$kitlg$fitRs$init$modelII)

rg$kitlg$runs$init$modelIII <-
  hzar.doFit(rg$kitlg$fitRs$init$modelIII)

rg$kitlg$runs$init$modelIV <-
  hzar.doFit(rg$kitlg$fitRs$init$modelIV)

rg$kitlg$runs$init$modelV <-
  hzar.doFit(rg$kitlg$fitRs$init$modelV)

rg$kitlg$runs$init$modelVI <-
  hzar.doFit(rg$kitlg$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$kitlg$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$kitlg$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$kitlg$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$kitlg$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$kitlg$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$kitlg$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$kitlg$fitRs$chains <-
  lapply(rg$kitlg$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$kitlg$fitRs$chains <-
  hzar.multiFitRequest(rg$kitlg$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$kitlg$runs$chains <-  hzar.doChain.multi(rg$kitlg$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$kitlg$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$kitlg$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$kitlg$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$kitlg$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$kitlg$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$kitlg$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$kitlg$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$kitlg$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$kitlg$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$kitlg$runs$init$modelI)
rg$kitlg$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$kitlg$runs$init$modelII)
rg$kitlg$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$kitlg$runs$init$modelIII)
rg$kitlg$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$kitlg$runs$init$modelIV)
rg$kitlg$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$kitlg$runs$init$modelV)
rg$kitlg$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$kitlg$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$kitlg$analysis$oDG <-
  hzar.make.obsDataGroup(rg$kitlg$analysis$initDGs)
rg$kitlg$analysis$oDG <-
  hzar.copyModelLabels(rg$kitlg$analysis$initDGs,
                       rg$kitlg$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$kitlg$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$kitlg$runs$chains,
                                hzar.dataGroup.add),
                         rg$kitlg$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$kitlg$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$kitlg$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$kitlg$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$kitlg$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$kitlg$analysis$model.name <-
        rownames(rg$kitlg$analysis$AICcTable
        )[[ which.min(rg$kitlg$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$kitlg$analysis$model.selected <-
  rg$kitlg$analysis$oDG$data.groups[[rg$kitlg$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$kitlg$analysis$model.selected,
                         names(rg$kitlg$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$kitlg$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$kitlg$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$kitlg$analysis$model.selected);

hzar.plot.fzCline(rg$back$analysis$model.selected);
hzar.plot.fzCline(rg$kitlg$analysis$model.selected);
## End Analysis

#dev.off()

### Analysis on PLXNC1 region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$plxnc1 <- list();
## Space to hold the observed data
rg$plxnc1$obs <- list();
## Space to hold the models to fit
rg$plxnc1$models <- list();
## Space to hold the compiled fit requests
rg$plxnc1$fitRs <- list();
## Space to hold the output data chains
rg$plxnc1$runs <- list();
## Space to hold the analysed data
rg$plxnc1$analysis <- list();

## Assign data for plxnc1 region
rg$plxnc1$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_plxnc1_mean,
                             rg.dists$hi_plxnc1_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$plxnc1$obs);

## Make a helper function to set cline models
rg.loadplxnc1model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$plxnc1$models[[id]] <<- hzar.makeCline1DFreq(rg$plxnc1$obs, scaling, tails)

rg.loadplxnc1model("fixed","none","modelI");
rg.loadplxnc1model("free" ,"none","modelII");
rg.loadplxnc1model("free" ,"both","modelIII");
rg.loadplxnc1model("free" ,"right","modelIV");
rg.loadplxnc1model("free" ,"left","modelV");
rg.loadplxnc1model("free" ,"mirror","modelVI");

## Check the default settings
print(rg$plxnc1$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$plxnc1$models <- sapply(rg$plxnc1$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$plxnc1$models)

## Compile each of the models to prepare for fitting
rg$plxnc1$fitRs$init <- sapply(rg$plxnc1$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$plxnc1$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$plxnc1$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$plxnc1$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$plxnc1$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$plxnc1$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$plxnc1$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$plxnc1$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$plxnc1$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$plxnc1$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$plxnc1$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$plxnc1$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$plxnc1$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$plxnc1$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$plxnc1$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$plxnc1$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$plxnc1$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$plxnc1$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$plxnc1$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$plxnc1$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$plxnc1$fitRs$init)

## Run the models for an initial chain
rg$plxnc1$runs$init <- list()

rg$plxnc1$runs$init$modelI <-
  hzar.doFit(rg$plxnc1$fitRs$init$modelI)

rg$plxnc1$runs$init$modelII <-
  hzar.doFit(rg$plxnc1$fitRs$init$modelII)

rg$plxnc1$runs$init$modelIII <-
  hzar.doFit(rg$plxnc1$fitRs$init$modelIII)

rg$plxnc1$runs$init$modelIV <-
  hzar.doFit(rg$plxnc1$fitRs$init$modelIV)

rg$plxnc1$runs$init$modelV <-
  hzar.doFit(rg$plxnc1$fitRs$init$modelV)

rg$plxnc1$runs$init$modelVI <-
  hzar.doFit(rg$plxnc1$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$plxnc1$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$plxnc1$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$plxnc1$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$plxnc1$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$plxnc1$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$plxnc1$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$plxnc1$fitRs$chains <-
  lapply(rg$plxnc1$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$plxnc1$fitRs$chains <-
  hzar.multiFitRequest(rg$plxnc1$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$plxnc1$runs$chains <-  hzar.doChain.multi(rg$plxnc1$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$plxnc1$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$plxnc1$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$plxnc1$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$plxnc1$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$plxnc1$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$plxnc1$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$plxnc1$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$plxnc1$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$plxnc1$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$plxnc1$runs$init$modelI)
rg$plxnc1$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$plxnc1$runs$init$modelII)
rg$plxnc1$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$plxnc1$runs$init$modelIII)
rg$plxnc1$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$plxnc1$runs$init$modelIV)
rg$plxnc1$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$plxnc1$runs$init$modelV)
rg$plxnc1$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$plxnc1$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$plxnc1$analysis$oDG <-
  hzar.make.obsDataGroup(rg$plxnc1$analysis$initDGs)
rg$plxnc1$analysis$oDG <-
  hzar.copyModelLabels(rg$plxnc1$analysis$initDGs,
                       rg$plxnc1$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$plxnc1$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$plxnc1$runs$chains,
                                hzar.dataGroup.add),
                         rg$plxnc1$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$plxnc1$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$plxnc1$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$plxnc1$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$plxnc1$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$plxnc1$analysis$model.name <-
        rownames(rg$plxnc1$analysis$AICcTable
        )[[ which.min(rg$plxnc1$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$plxnc1$analysis$model.selected <-
  rg$plxnc1$analysis$oDG$data.groups[[rg$plxnc1$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$plxnc1$analysis$model.selected,
                         names(rg$plxnc1$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$plxnc1$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$plxnc1$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$plxnc1$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on SPEF2-PRLR region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$spef2.prlr <- list();
## Space to hold the observed data
rg$spef2.prlr$obs <- list();
## Space to hold the models to fit
rg$spef2.prlr$models <- list();
## Space to hold the compiled fit requests
rg$spef2.prlr$fitRs <- list();
## Space to hold the output data chains
rg$spef2.prlr$runs <- list();
## Space to hold the analysed data
rg$spef2.prlr$analysis <- list();

## Assign data for spef2.prlr region
rg$spef2.prlr$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_spef2_prlr_mean,
                             rg.dists$hi_spef2_prlr_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$spef2.prlr$obs);

## Make a helper function to set cline models
rg.loadspef2.prlrmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$spef2.prlr$models[[id]] <<- hzar.makeCline1DFreq(rg$spef2.prlr$obs, scaling, tails)

rg.loadspef2.prlrmodel("fixed","none","modelI");
rg.loadspef2.prlrmodel("free" ,"none","modelII");
rg.loadspef2.prlrmodel("free" ,"both","modelIII");
rg.loadspef2.prlrmodel("free" ,"right","modelIV");
rg.loadspef2.prlrmodel("free" ,"left","modelV");
rg.loadspef2.prlrmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rg$spef2.prlr$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$spef2.prlr$models <- sapply(rg$spef2.prlr$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$spef2.prlr$models)

## Compile each of the models to prepare for fitting
rg$spef2.prlr$fitRs$init <- sapply(rg$spef2.prlr$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$spef2.prlr$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$spef2.prlr$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$spef2.prlr$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$spef2.prlr$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$spef2.prlr$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$spef2.prlr$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$spef2.prlr$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$spef2.prlr$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$spef2.prlr$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$spef2.prlr$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$spef2.prlr$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$spef2.prlr$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$spef2.prlr$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$spef2.prlr$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$spef2.prlr$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$spef2.prlr$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$spef2.prlr$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$spef2.prlr$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$spef2.prlr$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$spef2.prlr$fitRs$init)

## Run the models for an initial chain
rg$spef2.prlr$runs$init <- list()

rg$spef2.prlr$runs$init$modelI <-
  hzar.doFit(rg$spef2.prlr$fitRs$init$modelI)

rg$spef2.prlr$runs$init$modelII <-
  hzar.doFit(rg$spef2.prlr$fitRs$init$modelII)

rg$spef2.prlr$runs$init$modelIII <-
  hzar.doFit(rg$spef2.prlr$fitRs$init$modelIII)

rg$spef2.prlr$runs$init$modelIV <-
  hzar.doFit(rg$spef2.prlr$fitRs$init$modelIV)

rg$spef2.prlr$runs$init$modelV <-
  hzar.doFit(rg$spef2.prlr$fitRs$init$modelV)

rg$spef2.prlr$runs$init$modelVI <-
  hzar.doFit(rg$spef2.prlr$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$spef2.prlr$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$spef2.prlr$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$spef2.prlr$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$spef2.prlr$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$spef2.prlr$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$spef2.prlr$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$spef2.prlr$fitRs$chains <-
  lapply(rg$spef2.prlr$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$spef2.prlr$fitRs$chains <-
  hzar.multiFitRequest(rg$spef2.prlr$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$spef2.prlr$runs$chains <-  hzar.doChain.multi(rg$spef2.prlr$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$spef2.prlr$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$spef2.prlr$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$spef2.prlr$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$spef2.prlr$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$spef2.prlr$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$spef2.prlr$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$spef2.prlr$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$spef2.prlr$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$spef2.prlr$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$spef2.prlr$runs$init$modelI)
rg$spef2.prlr$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$spef2.prlr$runs$init$modelII)
rg$spef2.prlr$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$spef2.prlr$runs$init$modelIII)
rg$spef2.prlr$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$spef2.prlr$runs$init$modelIV)
rg$spef2.prlr$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$spef2.prlr$runs$init$modelV)
rg$spef2.prlr$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$spef2.prlr$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$spef2.prlr$analysis$oDG <-
  hzar.make.obsDataGroup(rg$spef2.prlr$analysis$initDGs)
rg$spef2.prlr$analysis$oDG <-
  hzar.copyModelLabels(rg$spef2.prlr$analysis$initDGs,
                       rg$spef2.prlr$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$spef2.prlr$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$spef2.prlr$runs$chains,
                                hzar.dataGroup.add),
                         rg$spef2.prlr$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$spef2.prlr$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$spef2.prlr$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$spef2.prlr$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$spef2.prlr$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$spef2.prlr$analysis$model.name <-
        rownames(rg$spef2.prlr$analysis$AICcTable
        )[[ which.min(rg$spef2.prlr$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$spef2.prlr$analysis$model.selected <-
  rg$spef2.prlr$analysis$oDG$data.groups[[rg$spef2.prlr$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$spef2.prlr$analysis$model.selected,
                         names(rg$spef2.prlr$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$spef2.prlr$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$spef2.prlr$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$spef2.prlr$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on SLC45A2 region (rustica-gutturalis)--------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$slc45a2 <- list();
## Space to hold the observed data
rg$slc45a2$obs <- list();
## Space to hold the models to fit
rg$slc45a2$models <- list();
## Space to hold the compiled fit requests
rg$slc45a2$fitRs <- list();
## Space to hold the output data chains
rg$slc45a2$runs <- list();
## Space to hold the analysed data
rg$slc45a2$analysis <- list();

## Assign data for slc45a2 region
rg$slc45a2$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_slc45a2_mean,
                             rg.dists$hi_slc45a2_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$slc45a2$obs);

## Make a helper function to set cline models
rg.loadslc45a2model <- function(scaling,tails,
                                id=paste(scaling,tails,sep="."))
  rg$slc45a2$models[[id]] <<- hzar.makeCline1DFreq(rg$slc45a2$obs, scaling, tails)

rg.loadslc45a2model("fixed","none","modelI");
rg.loadslc45a2model("free" ,"none","modelII");
rg.loadslc45a2model("free" ,"both","modelIII");
rg.loadslc45a2model("free" ,"right","modelIV");
rg.loadslc45a2model("free" ,"left","modelV");
rg.loadslc45a2model("free" ,"mirror","modelVI");

## Check the default settings
print(rg$slc45a2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$slc45a2$models <- sapply(rg$slc45a2$models,
                            hzar.model.addBoxReq,
                            -30 , 3526,
                            simplify=FALSE)

## Check the updated settings
print(rg$slc45a2$models)

## Compile each of the models to prepare for fitting
rg$slc45a2$fitRs$init <- sapply(rg$slc45a2$models,
                                hzar.first.fitRequest.old.ML,
                                obsData=rg$slc45a2$obs,
                                verbose=FALSE,
                                simplify=FALSE)

## Update the settings for the fitter if desired.
rg$slc45a2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$slc45a2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$slc45a2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$slc45a2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$slc45a2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$slc45a2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$slc45a2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$slc45a2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$slc45a2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$slc45a2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$slc45a2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$slc45a2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$slc45a2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$slc45a2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$slc45a2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$slc45a2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$slc45a2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$slc45a2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$slc45a2$fitRs$init)

## Run the models for an initial chain
rg$slc45a2$runs$init <- list()

rg$slc45a2$runs$init$modelI <-
  hzar.doFit(rg$slc45a2$fitRs$init$modelI)

rg$slc45a2$runs$init$modelII <-
  hzar.doFit(rg$slc45a2$fitRs$init$modelII)

rg$slc45a2$runs$init$modelIII <-
  hzar.doFit(rg$slc45a2$fitRs$init$modelIII)

rg$slc45a2$runs$init$modelIV <-
  hzar.doFit(rg$slc45a2$fitRs$init$modelIV)

rg$slc45a2$runs$init$modelV <-
  hzar.doFit(rg$slc45a2$fitRs$init$modelV)

rg$slc45a2$runs$init$modelVI <-
  hzar.doFit(rg$slc45a2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$slc45a2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$slc45a2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$slc45a2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$slc45a2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$slc45a2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$slc45a2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$slc45a2$fitRs$chains <-
  lapply(rg$slc45a2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$slc45a2$fitRs$chains <-
  hzar.multiFitRequest(rg$slc45a2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$slc45a2$runs$chains <-  hzar.doChain.multi(rg$slc45a2$fitRs$chains,
                                              doPar=TRUE,
                                              inOrder=FALSE,
                                              count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$slc45a2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$slc45a2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$slc45a2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$slc45a2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$slc45a2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$slc45a2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$slc45a2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$slc45a2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$slc45a2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$slc45a2$runs$init$modelI)
rg$slc45a2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$slc45a2$runs$init$modelII)
rg$slc45a2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$slc45a2$runs$init$modelIII)
rg$slc45a2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$slc45a2$runs$init$modelIV)
rg$slc45a2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$slc45a2$runs$init$modelV)
rg$slc45a2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$slc45a2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$slc45a2$analysis$oDG <-
  hzar.make.obsDataGroup(rg$slc45a2$analysis$initDGs)
rg$slc45a2$analysis$oDG <-
  hzar.copyModelLabels(rg$slc45a2$analysis$initDGs,
                       rg$slc45a2$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$slc45a2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$slc45a2$runs$chains,
                                hzar.dataGroup.add),
                         rg$slc45a2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$slc45a2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$slc45a2$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$slc45a2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$slc45a2$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$slc45a2$analysis$model.name <-
        rownames(rg$slc45a2$analysis$AICcTable
        )[[ which.min(rg$slc45a2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$slc45a2$analysis$model.selected <-
  rg$slc45a2$analysis$oDG$data.groups[[rg$slc45a2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$slc45a2$analysis$model.selected,
                         names(rg$slc45a2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$slc45a2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$slc45a2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$slc45a2$analysis$model.selected);
rect(xleft=rg$slc45a2$analysis$model.selected$ML.cline$param.all$center-((rg$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=rg$slc45a2$analysis$model.selected$ML.cline$param.all$center+((rg$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#8362AA')

## End Analysis

#dev.off()

### Analysis on BNC2 region (rustica-gutturalis)-----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$bnc2 <- list();
## Space to hold the observed data
rg$bnc2$obs <- list();
## Space to hold the models to fit
rg$bnc2$models <- list();
## Space to hold the compiled fit requests
rg$bnc2$fitRs <- list();
## Space to hold the output data chains
rg$bnc2$runs <- list();
## Space to hold the analysed data
rg$bnc2$analysis <- list();

## Assign data for bnc2 region
rg$bnc2$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_bnc2_mean,
                             rg.dists$hi_bnc2_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$bnc2$obs);

## Make a helper function to set cline models
rg.loadbnc2model <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  rg$bnc2$models[[id]] <<- hzar.makeCline1DFreq(rg$bnc2$obs, scaling, tails)

rg.loadbnc2model("fixed","none","modelI");
rg.loadbnc2model("free" ,"none","modelII");
rg.loadbnc2model("free" ,"both","modelIII");
rg.loadbnc2model("free" ,"right","modelIV");
rg.loadbnc2model("free" ,"left","modelV");
rg.loadbnc2model("free" ,"mirror","modelVI");

## Check the default settings
print(rg$bnc2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$bnc2$models <- sapply(rg$bnc2$models,
                         hzar.model.addBoxReq,
                         -30 , 3526,
                         simplify=FALSE)

## Check the updated settings
print(rg$bnc2$models)

## Compile each of the models to prepare for fitting
rg$bnc2$fitRs$init <- sapply(rg$bnc2$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=rg$bnc2$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
rg$bnc2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$bnc2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$bnc2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$bnc2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$bnc2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$bnc2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$bnc2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$bnc2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$bnc2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$bnc2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$bnc2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$bnc2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$bnc2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$bnc2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$bnc2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$bnc2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$bnc2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$bnc2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$bnc2$fitRs$init)

## Run the models for an initial chain
rg$bnc2$runs$init <- list()

rg$bnc2$runs$init$modelI <-
  hzar.doFit(rg$bnc2$fitRs$init$modelI)

rg$bnc2$runs$init$modelII <-
  hzar.doFit(rg$bnc2$fitRs$init$modelII)

rg$bnc2$runs$init$modelIII <-
  hzar.doFit(rg$bnc2$fitRs$init$modelIII)

rg$bnc2$runs$init$modelIV <-
  hzar.doFit(rg$bnc2$fitRs$init$modelIV)

rg$bnc2$runs$init$modelV <-
  hzar.doFit(rg$bnc2$fitRs$init$modelV)

rg$bnc2$runs$init$modelVI <-
  hzar.doFit(rg$bnc2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$bnc2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$bnc2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$bnc2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$bnc2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$bnc2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$bnc2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$bnc2$fitRs$chains <-
  lapply(rg$bnc2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$bnc2$fitRs$chains <-
  hzar.multiFitRequest(rg$bnc2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$bnc2$runs$chains <-  hzar.doChain.multi(rg$bnc2$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$bnc2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$bnc2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$bnc2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$bnc2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$bnc2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$bnc2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$bnc2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$bnc2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$bnc2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$bnc2$runs$init$modelI)
rg$bnc2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$bnc2$runs$init$modelII)
rg$bnc2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$bnc2$runs$init$modelIII)
rg$bnc2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$bnc2$runs$init$modelIV)
rg$bnc2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$bnc2$runs$init$modelV)
rg$bnc2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$bnc2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$bnc2$analysis$oDG <-
  hzar.make.obsDataGroup(rg$bnc2$analysis$initDGs)
rg$bnc2$analysis$oDG <-
  hzar.copyModelLabels(rg$bnc2$analysis$initDGs,
                       rg$bnc2$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$bnc2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$bnc2$runs$chains,
                                hzar.dataGroup.add),
                         rg$bnc2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$bnc2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$bnc2$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$bnc2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$bnc2$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$bnc2$analysis$model.name <-
        rownames(rg$bnc2$analysis$AICcTable
        )[[ which.min(rg$bnc2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$bnc2$analysis$model.selected <-
  rg$bnc2$analysis$oDG$data.groups[[rg$bnc2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$bnc2$analysis$model.selected,
                         names(rg$bnc2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$bnc2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$bnc2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$bnc2$analysis$model.selected);
rect(xleft=rg$bnc2$analysis$model.selected$ML.cline$param.all$center-((rg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=rg$bnc2$analysis$model.selected$ML.cline$param.all$center+((rg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#8362AA')

## End Analysis

#dev.off()

### Analysis on GNAQ region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$gnaq <- list();
## Space to hold the observed data
rg$gnaq$obs <- list();
## Space to hold the models to fit
rg$gnaq$models <- list();
## Space to hold the compiled fit requests
rg$gnaq$fitRs <- list();
## Space to hold the output data chains
rg$gnaq$runs <- list();
## Space to hold the analysed data
rg$gnaq$analysis <- list();

## Assign data for gnaq region
rg$gnaq$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_gnaq_mean,
                             rg.dists$hi_gnaq_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$gnaq$obs);

## Make a helper function to set cline models
rg.loadgnaqmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$gnaq$models[[id]] <<- hzar.makeCline1DFreq(rg$gnaq$obs, scaling, tails)

rg.loadgnaqmodel("fixed","none","modelI");
rg.loadgnaqmodel("free" ,"none","modelII");
rg.loadgnaqmodel("free" ,"both","modelIII");
rg.loadgnaqmodel("free" ,"right","modelIV");
rg.loadgnaqmodel("free" ,"left","modelV");
rg.loadgnaqmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rg$gnaq$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$gnaq$models <- sapply(rg$gnaq$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$gnaq$models)

## Compile each of the models to prepare for fitting
rg$gnaq$fitRs$init <- sapply(rg$gnaq$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$gnaq$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$gnaq$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$gnaq$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$gnaq$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$gnaq$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$gnaq$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$gnaq$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$gnaq$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$gnaq$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$gnaq$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$gnaq$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$gnaq$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$gnaq$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$gnaq$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$gnaq$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$gnaq$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$gnaq$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$gnaq$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$gnaq$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$gnaq$fitRs$init)

## Run the models for an initial chain
rg$gnaq$runs$init <- list()

rg$gnaq$runs$init$modelI <-
  hzar.doFit(rg$gnaq$fitRs$init$modelI)

rg$gnaq$runs$init$modelII <-
  hzar.doFit(rg$gnaq$fitRs$init$modelII)

rg$gnaq$runs$init$modelIII <-
  hzar.doFit(rg$gnaq$fitRs$init$modelIII)

rg$gnaq$runs$init$modelIV <-
  hzar.doFit(rg$gnaq$fitRs$init$modelIV)

rg$gnaq$runs$init$modelV <-
  hzar.doFit(rg$gnaq$fitRs$init$modelV)

rg$gnaq$runs$init$modelVI <-
  hzar.doFit(rg$gnaq$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$gnaq$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$gnaq$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$gnaq$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$gnaq$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$gnaq$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$gnaq$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$gnaq$fitRs$chains <-
  lapply(rg$gnaq$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$gnaq$fitRs$chains <-
  hzar.multiFitRequest(rg$gnaq$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$gnaq$runs$chains <-  hzar.doChain.multi(rg$gnaq$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$gnaq$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$gnaq$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$gnaq$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$gnaq$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$gnaq$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$gnaq$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$gnaq$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$gnaq$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$gnaq$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$gnaq$runs$init$modelI)
rg$gnaq$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$gnaq$runs$init$modelII)
rg$gnaq$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$gnaq$runs$init$modelIII)
rg$gnaq$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$gnaq$runs$init$modelIV)
rg$gnaq$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$gnaq$runs$init$modelV)
rg$gnaq$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$gnaq$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$gnaq$analysis$oDG <-
  hzar.make.obsDataGroup(rg$gnaq$analysis$initDGs)
rg$gnaq$analysis$oDG <-
  hzar.copyModelLabels(rg$gnaq$analysis$initDGs,
                       rg$gnaq$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$gnaq$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$gnaq$runs$chains,
                                hzar.dataGroup.add),
                         rg$gnaq$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$gnaq$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$gnaq$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$gnaq$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$gnaq$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$gnaq$analysis$model.name <-
        rownames(rg$gnaq$analysis$AICcTable
        )[[ which.min(rg$gnaq$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$gnaq$analysis$model.selected <-
  rg$gnaq$analysis$oDG$data.groups[[rg$gnaq$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$gnaq$analysis$model.selected,
                         names(rg$gnaq$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$gnaq$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$gnaq$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$gnaq$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on ROR2 region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$ror2 <- list();
## Space to hold the observed data
rg$ror2$obs <- list();
## Space to hold the models to fit
rg$ror2$models <- list();
## Space to hold the compiled fit requests
rg$ror2$fitRs <- list();
## Space to hold the output data chains
rg$ror2$runs <- list();
## Space to hold the analysed data
rg$ror2$analysis <- list();

## Assign data for ror2 region
rg$ror2$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_ror2_mean,
                             rg.dists$hi_ror2_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$ror2$obs);

## Make a helper function to set cline models
rg.loadror2model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$ror2$models[[id]] <<- hzar.makeCline1DFreq(rg$ror2$obs, scaling, tails)

rg.loadror2model("fixed","none","modelI");
rg.loadror2model("free" ,"none","modelII");
rg.loadror2model("free" ,"both","modelIII");
rg.loadror2model("free" ,"right","modelIV");
rg.loadror2model("free" ,"left","modelV");
rg.loadror2model("free" ,"mirror","modelVI");

## Check the default settings
print(rg$ror2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$ror2$models <- sapply(rg$ror2$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$ror2$models)

## Compile each of the models to prepare for fitting
rg$ror2$fitRs$init <- sapply(rg$ror2$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$ror2$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$ror2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ror2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ror2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$ror2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ror2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ror2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$ror2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ror2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ror2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$ror2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ror2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ror2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$ror2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ror2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ror2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$ror2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ror2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ror2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$ror2$fitRs$init)

## Run the models for an initial chain
rg$ror2$runs$init <- list()

rg$ror2$runs$init$modelI <-
  hzar.doFit(rg$ror2$fitRs$init$modelI)

rg$ror2$runs$init$modelII <-
  hzar.doFit(rg$ror2$fitRs$init$modelII)

rg$ror2$runs$init$modelIII <-
  hzar.doFit(rg$ror2$fitRs$init$modelIII)

rg$ror2$runs$init$modelIV <-
  hzar.doFit(rg$ror2$fitRs$init$modelIV)

rg$ror2$runs$init$modelV <-
  hzar.doFit(rg$ror2$fitRs$init$modelV)

rg$ror2$runs$init$modelVI <-
  hzar.doFit(rg$ror2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$ror2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$ror2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$ror2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$ror2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$ror2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$ror2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$ror2$fitRs$chains <-
  lapply(rg$ror2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$ror2$fitRs$chains <-
  hzar.multiFitRequest(rg$ror2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$ror2$runs$chains <-  hzar.doChain.multi(rg$ror2$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$ror2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$ror2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$ror2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$ror2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$ror2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$ror2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$ror2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$ror2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$ror2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$ror2$runs$init$modelI)
rg$ror2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$ror2$runs$init$modelII)
rg$ror2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$ror2$runs$init$modelIII)
rg$ror2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$ror2$runs$init$modelIV)
rg$ror2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$ror2$runs$init$modelV)
rg$ror2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$ror2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$ror2$analysis$oDG <-
  hzar.make.obsDataGroup(rg$ror2$analysis$initDGs)
rg$ror2$analysis$oDG <-
  hzar.copyModelLabels(rg$ror2$analysis$initDGs,
                       rg$ror2$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$ror2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$ror2$runs$chains,
                                hzar.dataGroup.add),
                         rg$ror2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$ror2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$ror2$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$ror2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$ror2$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$ror2$analysis$model.name <-
        rownames(rg$ror2$analysis$AICcTable
        )[[ which.min(rg$ror2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$ror2$analysis$model.selected <-
  rg$ror2$analysis$oDG$data.groups[[rg$ror2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$ror2$analysis$model.selected,
                         names(rg$ror2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$ror2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$ror2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$ror2$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on APC-CAMK4 region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$apc.camk4 <- list();
## Space to hold the observed data
rg$apc.camk4$obs <- list();
## Space to hold the models to fit
rg$apc.camk4$models <- list();
## Space to hold the compiled fit requests
rg$apc.camk4$fitRs <- list();
## Space to hold the output data chains
rg$apc.camk4$runs <- list();
## Space to hold the analysed data
rg$apc.camk4$analysis <- list();

## Assign data for apc.camk4 region
rg$apc.camk4$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_apc_camk4_mean,
                             rg.dists$hi_apc_camk4_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$apc.camk4$obs);

## Make a helper function to set cline models
rg.loadapc.camk4model <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$apc.camk4$models[[id]] <<- hzar.makeCline1DFreq(rg$apc.camk4$obs, scaling, tails)

rg.loadapc.camk4model("fixed","none","modelI");
rg.loadapc.camk4model("free" ,"none","modelII");
rg.loadapc.camk4model("free" ,"both","modelIII");
rg.loadapc.camk4model("free" ,"right","modelIV");
rg.loadapc.camk4model("free" ,"left","modelV");
rg.loadapc.camk4model("free" ,"mirror","modelVI");

## Check the default settings
print(rg$apc.camk4$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$apc.camk4$models <- sapply(rg$apc.camk4$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$apc.camk4$models)

## Compile each of the models to prepare for fitting
rg$apc.camk4$fitRs$init <- sapply(rg$apc.camk4$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$apc.camk4$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$apc.camk4$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$apc.camk4$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$apc.camk4$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$apc.camk4$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$apc.camk4$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$apc.camk4$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$apc.camk4$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$apc.camk4$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$apc.camk4$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$apc.camk4$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$apc.camk4$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$apc.camk4$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$apc.camk4$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$apc.camk4$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$apc.camk4$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$apc.camk4$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$apc.camk4$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$apc.camk4$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$apc.camk4$fitRs$init)

## Run the models for an initial chain
rg$apc.camk4$runs$init <- list()

rg$apc.camk4$runs$init$modelI <-
  hzar.doFit(rg$apc.camk4$fitRs$init$modelI)

rg$apc.camk4$runs$init$modelII <-
  hzar.doFit(rg$apc.camk4$fitRs$init$modelII)

rg$apc.camk4$runs$init$modelIII <-
  hzar.doFit(rg$apc.camk4$fitRs$init$modelIII)

rg$apc.camk4$runs$init$modelIV <-
  hzar.doFit(rg$apc.camk4$fitRs$init$modelIV)

rg$apc.camk4$runs$init$modelV <-
  hzar.doFit(rg$apc.camk4$fitRs$init$modelV)

rg$apc.camk4$runs$init$modelVI <-
  hzar.doFit(rg$apc.camk4$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$apc.camk4$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$apc.camk4$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$apc.camk4$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$apc.camk4$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$apc.camk4$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$apc.camk4$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$apc.camk4$fitRs$chains <-
  lapply(rg$apc.camk4$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$apc.camk4$fitRs$chains <-
  hzar.multiFitRequest(rg$apc.camk4$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$apc.camk4$runs$chains <-  hzar.doChain.multi(rg$apc.camk4$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$apc.camk4$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$apc.camk4$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$apc.camk4$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$apc.camk4$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$apc.camk4$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$apc.camk4$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$apc.camk4$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$apc.camk4$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$apc.camk4$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$apc.camk4$runs$init$modelI)
rg$apc.camk4$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$apc.camk4$runs$init$modelII)
rg$apc.camk4$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$apc.camk4$runs$init$modelIII)
rg$apc.camk4$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$apc.camk4$runs$init$modelIV)
rg$apc.camk4$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$apc.camk4$runs$init$modelV)
rg$apc.camk4$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$apc.camk4$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$apc.camk4$analysis$oDG <-
  hzar.make.obsDataGroup(rg$apc.camk4$analysis$initDGs)
rg$apc.camk4$analysis$oDG <-
  hzar.copyModelLabels(rg$apc.camk4$analysis$initDGs,
                       rg$apc.camk4$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$apc.camk4$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$apc.camk4$runs$chains,
                                hzar.dataGroup.add),
                         rg$apc.camk4$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$apc.camk4$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$apc.camk4$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$apc.camk4$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$apc.camk4$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$apc.camk4$analysis$model.name <-
        rownames(rg$apc.camk4$analysis$AICcTable
        )[[ which.min(rg$apc.camk4$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$apc.camk4$analysis$model.selected <-
  rg$apc.camk4$analysis$oDG$data.groups[[rg$apc.camk4$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$apc.camk4$analysis$model.selected,
                         names(rg$apc.camk4$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$apc.camk4$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$apc.camk4$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$apc.camk4$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on ICE1-lncRNA region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$ice1.lncrna <- list();
## Space to hold the observed data
rg$ice1.lncrna$obs <- list();
## Space to hold the models to fit
rg$ice1.lncrna$models <- list();
## Space to hold the compiled fit requests
rg$ice1.lncrna$fitRs <- list();
## Space to hold the output data chains
rg$ice1.lncrna$runs <- list();
## Space to hold the analysed data
rg$ice1.lncrna$analysis <- list();

## Assign data for ice1.lncrna region
rg$ice1.lncrna$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_ice1_lncrna_mean,
                             rg.dists$hi_ice1_lncrna_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$ice1.lncrna$obs);

## Make a helper function to set cline models
rg.loadice1.lncrnamodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$ice1.lncrna$models[[id]] <<- hzar.makeCline1DFreq(rg$ice1.lncrna$obs, scaling, tails)

rg.loadice1.lncrnamodel("fixed","none","modelI");
rg.loadice1.lncrnamodel("free" ,"none","modelII");
rg.loadice1.lncrnamodel("free" ,"both","modelIII");
rg.loadice1.lncrnamodel("free" ,"right","modelIV");
rg.loadice1.lncrnamodel("free" ,"left","modelV");
rg.loadice1.lncrnamodel("free" ,"mirror","modelVI");

## Check the default settings
print(rg$ice1.lncrna$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$ice1.lncrna$models <- sapply(rg$ice1.lncrna$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$ice1.lncrna$models)

## Compile each of the models to prepare for fitting
rg$ice1.lncrna$fitRs$init <- sapply(rg$ice1.lncrna$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$ice1.lncrna$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$ice1.lncrna$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ice1.lncrna$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ice1.lncrna$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$ice1.lncrna$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ice1.lncrna$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ice1.lncrna$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$ice1.lncrna$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ice1.lncrna$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ice1.lncrna$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$ice1.lncrna$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ice1.lncrna$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ice1.lncrna$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$ice1.lncrna$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ice1.lncrna$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ice1.lncrna$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$ice1.lncrna$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$ice1.lncrna$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$ice1.lncrna$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$ice1.lncrna$fitRs$init)

## Run the models for an initial chain
rg$ice1.lncrna$runs$init <- list()

rg$ice1.lncrna$runs$init$modelI <-
  hzar.doFit(rg$ice1.lncrna$fitRs$init$modelI)

rg$ice1.lncrna$runs$init$modelII <-
  hzar.doFit(rg$ice1.lncrna$fitRs$init$modelII)

rg$ice1.lncrna$runs$init$modelIII <-
  hzar.doFit(rg$ice1.lncrna$fitRs$init$modelIII)

rg$ice1.lncrna$runs$init$modelIV <-
  hzar.doFit(rg$ice1.lncrna$fitRs$init$modelIV)

rg$ice1.lncrna$runs$init$modelV <-
  hzar.doFit(rg$ice1.lncrna$fitRs$init$modelV)

rg$ice1.lncrna$runs$init$modelVI <-
  hzar.doFit(rg$ice1.lncrna$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$ice1.lncrna$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$ice1.lncrna$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$ice1.lncrna$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$ice1.lncrna$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$ice1.lncrna$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$ice1.lncrna$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$ice1.lncrna$fitRs$chains <-
  lapply(rg$ice1.lncrna$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$ice1.lncrna$fitRs$chains <-
  hzar.multiFitRequest(rg$ice1.lncrna$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$ice1.lncrna$runs$chains <-  hzar.doChain.multi(rg$ice1.lncrna$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$ice1.lncrna$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$ice1.lncrna$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$ice1.lncrna$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$ice1.lncrna$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$ice1.lncrna$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$ice1.lncrna$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$ice1.lncrna$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$ice1.lncrna$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$ice1.lncrna$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$ice1.lncrna$runs$init$modelI)
rg$ice1.lncrna$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$ice1.lncrna$runs$init$modelII)
rg$ice1.lncrna$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$ice1.lncrna$runs$init$modelIII)
rg$ice1.lncrna$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$ice1.lncrna$runs$init$modelIV)
rg$ice1.lncrna$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$ice1.lncrna$runs$init$modelV)
rg$ice1.lncrna$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$ice1.lncrna$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$ice1.lncrna$analysis$oDG <-
  hzar.make.obsDataGroup(rg$ice1.lncrna$analysis$initDGs)
rg$ice1.lncrna$analysis$oDG <-
  hzar.copyModelLabels(rg$ice1.lncrna$analysis$initDGs,
                       rg$ice1.lncrna$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$ice1.lncrna$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$ice1.lncrna$runs$chains,
                                hzar.dataGroup.add),
                         rg$ice1.lncrna$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$ice1.lncrna$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$ice1.lncrna$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$ice1.lncrna$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$ice1.lncrna$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$ice1.lncrna$analysis$model.name <-
        rownames(rg$ice1.lncrna$analysis$AICcTable
        )[[ which.min(rg$ice1.lncrna$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$ice1.lncrna$analysis$model.selected <-
  rg$ice1.lncrna$analysis$oDG$data.groups[[rg$ice1.lncrna$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$ice1.lncrna$analysis$model.selected,
                         names(rg$ice1.lncrna$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$ice1.lncrna$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$ice1.lncrna$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$ice1.lncrna$analysis$model.selected);

## End Analysis

#dev.off()


### Analysis on PDE1C region (rustica-gutturalis)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg$pde1c <- list();
## Space to hold the observed data
rg$pde1c$obs <- list();
## Space to hold the models to fit
rg$pde1c$models <- list();
## Space to hold the compiled fit requests
rg$pde1c$fitRs <- list();
## Space to hold the output data chains
rg$pde1c$runs <- list();
## Space to hold the analysed data
rg$pde1c$analysis <- list();

## Assign data for pde1c region
rg$pde1c$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$hi_pde1c_mean,
                             rg.dists$hi_pde1c_count);

## Look at a graph of the observed data
hzar.plot.obsData(rg$pde1c$obs);

## Make a helper function to set cline models
rg.loadpde1cmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rg$pde1c$models[[id]] <<- hzar.makeCline1DFreq(rg$pde1c$obs, scaling, tails)

rg.loadpde1cmodel("fixed","none","modelI");
rg.loadpde1cmodel("free" ,"none","modelII");
rg.loadpde1cmodel("free" ,"both","modelIII");
rg.loadpde1cmodel("free" ,"right","modelIV");
rg.loadpde1cmodel("free" ,"left","modelV");
rg.loadpde1cmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rg$pde1c$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 3495.754 km.

min(rg.dists$km)
max(rg.dists$km)

rg$pde1c$models <- sapply(rg$pde1c$models,
                          hzar.model.addBoxReq,
                          -30 , 3526,
                          simplify=FALSE)

## Check the updated settings
print(rg$pde1c$models)

## Compile each of the models to prepare for fitting
rg$pde1c$fitRs$init <- sapply(rg$pde1c$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rg$pde1c$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rg$pde1c$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$pde1c$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$pde1c$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg$pde1c$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$pde1c$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$pde1c$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg$pde1c$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$pde1c$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$pde1c$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$pde1c$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$pde1c$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$pde1c$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$pde1c$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$pde1c$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$pde1c$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg$pde1c$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg$pde1c$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg$pde1c$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg$pde1c$fitRs$init)

## Run the models for an initial chain
rg$pde1c$runs$init <- list()

rg$pde1c$runs$init$modelI <-
  hzar.doFit(rg$pde1c$fitRs$init$modelI)

rg$pde1c$runs$init$modelII <-
  hzar.doFit(rg$pde1c$fitRs$init$modelII)

rg$pde1c$runs$init$modelIII <-
  hzar.doFit(rg$pde1c$fitRs$init$modelIII)

rg$pde1c$runs$init$modelIV <-
  hzar.doFit(rg$pde1c$fitRs$init$modelIV)

rg$pde1c$runs$init$modelV <-
  hzar.doFit(rg$pde1c$fitRs$init$modelV)

rg$pde1c$runs$init$modelVI <-
  hzar.doFit(rg$pde1c$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rg$pde1c$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rg$pde1c$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rg$pde1c$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rg$pde1c$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rg$pde1c$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rg$pde1c$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rg$pde1c$fitRs$chains <-
  lapply(rg$pde1c$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg$pde1c$fitRs$chains <-
  hzar.multiFitRequest(rg$pde1c$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg$pde1c$runs$chains <-  hzar.doChain.multi(rg$pde1c$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg$pde1c$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg$pde1c$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg$pde1c$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg$pde1c$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg$pde1c$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg$pde1c$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg$pde1c$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg$pde1c$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg$pde1c$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg$pde1c$runs$init$modelI)
rg$pde1c$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg$pde1c$runs$init$modelII)
rg$pde1c$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg$pde1c$runs$init$modelIII)
rg$pde1c$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg$pde1c$runs$init$modelIV)
rg$pde1c$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg$pde1c$runs$init$modelV)
rg$pde1c$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg$pde1c$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg$pde1c$analysis$oDG <-
  hzar.make.obsDataGroup(rg$pde1c$analysis$initDGs)
rg$pde1c$analysis$oDG <-
  hzar.copyModelLabels(rg$pde1c$analysis$initDGs,
                       rg$pde1c$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg$pde1c$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg$pde1c$runs$chains,
                                hzar.dataGroup.add),
                         rg$pde1c$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg$pde1c$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rg$pde1c$analysis$oDG);

## Do model selection based on the AICc scores
print(rg$pde1c$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$pde1c$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg$pde1c$analysis$model.name <-
        rownames(rg$pde1c$analysis$AICcTable
        )[[ which.min(rg$pde1c$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg$pde1c$analysis$model.selected <-
  rg$pde1c$analysis$oDG$data.groups[[rg$pde1c$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$pde1c$analysis$model.selected,
                         names(rg$pde1c$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$pde1c$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$pde1c$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rg$pde1c$analysis$model.selected);

## End Analysis

#dev.off()

### Save analysis results to Rdata object-----------------------------------

save(rg, file = './candidate_results_hzar/hzar_data_rustica-gutturalis.RData')
