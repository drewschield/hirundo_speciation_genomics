############################################################################
# Barn swallow geographic clines in hybrid zone transects
############################################################################

# This script contains commands for analysis of individual hybrid index and
# geographic clines for 'outlier' regions of the genome versus background
# to understand if genomic regions associated with mate choice traits show
# restricted gene flow (i.e., steeper/narrower clines).

# This workflow involves:
# 1. Estimation of hybrid indices from outlier regions for the rustica-tytleri hybrid zone
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
load("./candidate_results_hzar/hzar_data_rustica-tytleri.RData")

## Introgress - estimation of hybrid index (rustica-tytleri)---------------

## Read in VCF data
rt_ru.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.rustica-tytleri_rustica.vcf.gz')
rt_ty.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.rustica-tytleri_tytleri.vcf.gz')
rt_rt.kitlg.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-KITLG.rustica-tytleri_hybrids.vcf.gz')

rt_ru.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.rustica-tytleri_rustica.vcf.gz')
rt_ty.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.rustica-tytleri_tytleri.vcf.gz')
rt_rt.plxnc1.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PLXNC1.rustica-tytleri_hybrids.vcf.gz')

rt_ru.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.rustica-tytleri_rustica.vcf.gz')
rt_ty.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.rustica-tytleri_tytleri.vcf.gz')
rt_rt.spef2.prlr.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SPEF2-PRLR.rustica-tytleri_hybrids.vcf.gz')

rt_ru.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.rustica-tytleri_rustica.vcf.gz')
rt_ty.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.rustica-tytleri_tytleri.vcf.gz')
rt_rt.slc45a2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-SLC45A2.rustica-tytleri_hybrids.vcf.gz')

rt_ru.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.rustica-tytleri_rustica.vcf.gz')
rt_ty.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.rustica-tytleri_tytleri.vcf.gz')
rt_rt.bnc2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-BNC2.rustica-tytleri_hybrids.vcf.gz')

rt_ru.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.rustica-tytleri_rustica.vcf.gz')
rt_ty.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.rustica-tytleri_tytleri.vcf.gz')
rt_rt.gnaq.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-GNAQ.rustica-tytleri_hybrids.vcf.gz')

rt_ru.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.rustica-tytleri_rustica.vcf.gz')
rt_ty.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.rustica-tytleri_tytleri.vcf.gz')
rt_rt.ror2.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ROR2.rustica-tytleri_hybrids.vcf.gz')

rt_ru.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.rustica-tytleri_rustica.vcf.gz')
rt_ty.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.rustica-tytleri_tytleri.vcf.gz')
rt_rt.apc.camk4.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-APC-CAMK4.rustica-tytleri_hybrids.vcf.gz')

rt_ru.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.rustica-tytleri_rustica.vcf.gz')
rt_ty.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.rustica-tytleri_tytleri.vcf.gz')
rt_rt.ice1.lncrna.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-ICE1-lncRNA.rustica-tytleri_hybrids.vcf.gz')

rt_ru.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.rustica-tytleri_rustica.vcf.gz')
rt_ty.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.rustica-tytleri_tytleri.vcf.gz')
rt_rt.pde1c.vcfr <- read.vcfR('./candidate_input_hzar/anc-info.snps.outlier-PDE1C.rustica-tytleri_hybrids.vcf.gz')

## Extract genotype matrices

# Candidate regions
rt_ru.kitlg.gt <- extract.gt(rt_ru.kitlg.vcfr)
rt_ty.kitlg.gt <- extract.gt(rt_ty.kitlg.vcfr)
rt_rt.kitlg.gt <- extract.gt(rt_rt.kitlg.vcfr)

rt_ru.plxnc1.gt <- extract.gt(rt_ru.plxnc1.vcfr)
rt_ty.plxnc1.gt <- extract.gt(rt_ty.plxnc1.vcfr)
rt_rt.plxnc1.gt <- extract.gt(rt_rt.plxnc1.vcfr)

rt_ru.spef2.prlr.gt <- extract.gt(rt_ru.spef2.prlr.vcfr)
rt_ty.spef2.prlr.gt <- extract.gt(rt_ty.spef2.prlr.vcfr)
rt_rt.spef2.prlr.gt <- extract.gt(rt_rt.spef2.prlr.vcfr)

rt_ru.slc45a2.gt <- extract.gt(rt_ru.slc45a2.vcfr)
rt_ty.slc45a2.gt <- extract.gt(rt_ty.slc45a2.vcfr)
rt_rt.slc45a2.gt <- extract.gt(rt_rt.slc45a2.vcfr)

rt_ru.bnc2.gt <- extract.gt(rt_ru.bnc2.vcfr)
rt_ty.bnc2.gt <- extract.gt(rt_ty.bnc2.vcfr)
rt_rt.bnc2.gt <- extract.gt(rt_rt.bnc2.vcfr)

rt_ru.gnaq.gt <- extract.gt(rt_ru.gnaq.vcfr)
rt_ty.gnaq.gt <- extract.gt(rt_ty.gnaq.vcfr)
rt_rt.gnaq.gt <- extract.gt(rt_rt.gnaq.vcfr)

rt_ru.ror2.gt <- extract.gt(rt_ru.ror2.vcfr)
rt_ty.ror2.gt <- extract.gt(rt_ty.ror2.vcfr)
rt_rt.ror2.gt <- extract.gt(rt_rt.ror2.vcfr)

rt_ru.apc.camk4.gt <- extract.gt(rt_ru.apc.camk4.vcfr)
rt_ty.apc.camk4.gt <- extract.gt(rt_ty.apc.camk4.vcfr)
rt_rt.apc.camk4.gt <- extract.gt(rt_rt.apc.camk4.vcfr)

rt_ru.ice1.lncrna.gt <- extract.gt(rt_ru.ice1.lncrna.vcfr)
rt_ty.ice1.lncrna.gt <- extract.gt(rt_ty.ice1.lncrna.vcfr)
rt_rt.ice1.lncrna.gt <- extract.gt(rt_rt.ice1.lncrna.vcfr)

rt_ru.pde1c.gt <- extract.gt(rt_ru.pde1c.vcfr)
rt_ty.pde1c.gt <- extract.gt(rt_ty.pde1c.vcfr)
rt_rt.pde1c.gt <- extract.gt(rt_rt.pde1c.vcfr)

## Set locus sets
rt.kitlg.loci.dat <- as.data.frame(rt_rt.kitlg.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.plxnc1.loci.dat <- as.data.frame(rt_rt.plxnc1.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.spef2.prlr.loci.dat <- as.data.frame(rt_rt.spef2.prlr.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.slc45a2.loci.dat <- as.data.frame(rt_rt.slc45a2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.bnc2.loci.dat <- as.data.frame(rt_rt.bnc2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.gnaq.loci.dat <- as.data.frame(rt_rt.gnaq.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.ror2.loci.dat <- as.data.frame(rt_rt.ror2.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.apc.camk4.loci.dat <- as.data.frame(rt_rt.apc.camk4.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.ice1.lncrna.loci.dat <- as.data.frame(rt_rt.ice1.lncrna.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

rt.pde1c.loci.dat <- as.data.frame(rt_rt.pde1c.gt) %>%
  rownames_to_column("locus") %>%
  dplyr::select(locus) %>%
  mutate(type = "c")

## Prepare data for analysis
rt_ru.kitlg.dat <- prepare.data(rt_ru.kitlg.gt, loci.data = rt.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.kitlg.gt, parental2 = rt_ty.kitlg.gt)
rt_ty.kitlg.dat <- prepare.data(rt_ty.kitlg.gt, loci.data = rt.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.kitlg.gt, parental2 = rt_ty.kitlg.gt)
rt_rt.kitlg.dat <- prepare.data(rt_rt.kitlg.gt, loci.data = rt.kitlg.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.kitlg.gt, parental2 = rt_ty.kitlg.gt)

rt_ru.plxnc1.dat <- prepare.data(rt_ru.plxnc1.gt, loci.data = rt.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.plxnc1.gt, parental2 = rt_ty.plxnc1.gt)
rt_ty.plxnc1.dat <- prepare.data(rt_ty.plxnc1.gt, loci.data = rt.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.plxnc1.gt, parental2 = rt_ty.plxnc1.gt)
rt_rt.plxnc1.dat <- prepare.data(rt_rt.plxnc1.gt, loci.data = rt.plxnc1.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.plxnc1.gt, parental2 = rt_ty.plxnc1.gt)

rt_ru.spef2.prlr.dat <- prepare.data(rt_ru.spef2.prlr.gt, loci.data = rt.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.spef2.prlr.gt, parental2 = rt_ty.spef2.prlr.gt)
rt_ty.spef2.prlr.dat <- prepare.data(rt_ty.spef2.prlr.gt, loci.data = rt.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.spef2.prlr.gt, parental2 = rt_ty.spef2.prlr.gt)
rt_rt.spef2.prlr.dat <- prepare.data(rt_rt.spef2.prlr.gt, loci.data = rt.spef2.prlr.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.spef2.prlr.gt, parental2 = rt_ty.spef2.prlr.gt)

rt_ru.slc45a2.dat <- prepare.data(rt_ru.slc45a2.gt, loci.data = rt.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.slc45a2.gt, parental2 = rt_ty.slc45a2.gt)
rt_ty.slc45a2.dat <- prepare.data(rt_ty.slc45a2.gt, loci.data = rt.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.slc45a2.gt, parental2 = rt_ty.slc45a2.gt)
rt_rt.slc45a2.dat <- prepare.data(rt_rt.slc45a2.gt, loci.data = rt.slc45a2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.slc45a2.gt, parental2 = rt_ty.slc45a2.gt)

rt_ru.bnc2.dat <- prepare.data(rt_ru.bnc2.gt, loci.data = rt.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.bnc2.gt, parental2 = rt_ty.bnc2.gt)
rt_ty.bnc2.dat <- prepare.data(rt_ty.bnc2.gt, loci.data = rt.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.bnc2.gt, parental2 = rt_ty.bnc2.gt)
rt_rt.bnc2.dat <- prepare.data(rt_rt.bnc2.gt, loci.data = rt.bnc2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.bnc2.gt, parental2 = rt_ty.bnc2.gt)

rt_ru.gnaq.dat <- prepare.data(rt_ru.gnaq.gt, loci.data = rt.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.gnaq.gt, parental2 = rt_ty.gnaq.gt)
rt_ty.gnaq.dat <- prepare.data(rt_ty.gnaq.gt, loci.data = rt.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.gnaq.gt, parental2 = rt_ty.gnaq.gt)
rt_rt.gnaq.dat <- prepare.data(rt_rt.gnaq.gt, loci.data = rt.gnaq.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.gnaq.gt, parental2 = rt_ty.gnaq.gt)

rt_ru.ror2.dat <- prepare.data(rt_ru.ror2.gt, loci.data = rt.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.ror2.gt, parental2 = rt_ty.ror2.gt)
rt_ty.ror2.dat <- prepare.data(rt_ty.ror2.gt, loci.data = rt.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.ror2.gt, parental2 = rt_ty.ror2.gt)
rt_rt.ror2.dat <- prepare.data(rt_rt.ror2.gt, loci.data = rt.ror2.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.ror2.gt, parental2 = rt_ty.ror2.gt)

rt_ru.apc.camk4.dat <- prepare.data(rt_ru.apc.camk4.gt, loci.data = rt.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.apc.camk4.gt, parental2 = rt_ty.apc.camk4.gt)
rt_ty.apc.camk4.dat <- prepare.data(rt_ty.apc.camk4.gt, loci.data = rt.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.apc.camk4.gt, parental2 = rt_ty.apc.camk4.gt)
rt_rt.apc.camk4.dat <- prepare.data(rt_rt.apc.camk4.gt, loci.data = rt.apc.camk4.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.apc.camk4.gt, parental2 = rt_ty.apc.camk4.gt)

rt_ru.ice1.lncrna.dat <- prepare.data(rt_ru.ice1.lncrna.gt, loci.data = rt.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.ice1.lncrna.gt, parental2 = rt_ty.ice1.lncrna.gt)
rt_ty.ice1.lncrna.dat <- prepare.data(rt_ty.ice1.lncrna.gt, loci.data = rt.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.ice1.lncrna.gt, parental2 = rt_ty.ice1.lncrna.gt)
rt_rt.ice1.lncrna.dat <- prepare.data(rt_rt.ice1.lncrna.gt, loci.data = rt.ice1.lncrna.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.ice1.lncrna.gt, parental2 = rt_ty.ice1.lncrna.gt)

rt_ru.pde1c.dat <- prepare.data(rt_ru.pde1c.gt, loci.data = rt.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.pde1c.gt, parental2 = rt_ty.pde1c.gt)
rt_ty.pde1c.dat <- prepare.data(rt_ty.pde1c.gt, loci.data = rt.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.pde1c.gt, parental2 = rt_ty.pde1c.gt)
rt_rt.pde1c.dat <- prepare.data(rt_rt.pde1c.gt, loci.data = rt.pde1c.loci.dat, fixed = FALSE, pop.id = FALSE, ind.id = FALSE, parental1 = rt_ru.pde1c.gt, parental2 = rt_ty.pde1c.gt)

## Estimate hybrid index
rt_ru.kitlg.h <- est.h(rt_ru.kitlg.dat, loci.data = rt.kitlg.loci.dat, fixed = FALSE)
rt_ty.kitlg.h <- est.h(rt_ty.kitlg.dat, loci.data = rt.kitlg.loci.dat, fixed = FALSE)
rt_rt.kitlg.h <- est.h(rt_rt.kitlg.dat, loci.data = rt.kitlg.loci.dat, fixed = FALSE)

rt_ru.plxnc1.h <- est.h(rt_ru.plxnc1.dat, loci.data = rt.plxnc1.loci.dat, fixed = FALSE)
rt_ty.plxnc1.h <- est.h(rt_ty.plxnc1.dat, loci.data = rt.plxnc1.loci.dat, fixed = FALSE)
rt_rt.plxnc1.h <- est.h(rt_rt.plxnc1.dat, loci.data = rt.plxnc1.loci.dat, fixed = FALSE)

rt_ru.spef2.prlr.h <- est.h(rt_ru.spef2.prlr.dat, loci.data = rt.spef2.prlr.loci.dat, fixed = FALSE)
rt_ty.spef2.prlr.h <- est.h(rt_ty.spef2.prlr.dat, loci.data = rt.spef2.prlr.loci.dat, fixed = FALSE)
rt_rt.spef2.prlr.h <- est.h(rt_rt.spef2.prlr.dat, loci.data = rt.spef2.prlr.loci.dat, fixed = FALSE)

rt_ru.slc45a2.h <- est.h(rt_ru.slc45a2.dat, loci.data = rt.slc45a2.loci.dat, fixed = FALSE)
rt_ty.slc45a2.h <- est.h(rt_ty.slc45a2.dat, loci.data = rt.slc45a2.loci.dat, fixed = FALSE)
rt_rt.slc45a2.h <- est.h(rt_rt.slc45a2.dat, loci.data = rt.slc45a2.loci.dat, fixed = FALSE)

rt_ru.bnc2.h <- est.h(rt_ru.bnc2.dat, loci.data = rt.bnc2.loci.dat, fixed = FALSE)
rt_ty.bnc2.h <- est.h(rt_ty.bnc2.dat, loci.data = rt.bnc2.loci.dat, fixed = FALSE)
rt_rt.bnc2.h <- est.h(rt_rt.bnc2.dat, loci.data = rt.bnc2.loci.dat, fixed = FALSE)

rt_ru.gnaq.h <- est.h(rt_ru.gnaq.dat, loci.data = rt.gnaq.loci.dat, fixed = FALSE)
rt_ty.gnaq.h <- est.h(rt_ty.gnaq.dat, loci.data = rt.gnaq.loci.dat, fixed = FALSE)
rt_rt.gnaq.h <- est.h(rt_rt.gnaq.dat, loci.data = rt.gnaq.loci.dat, fixed = FALSE)

rt_ru.ror2.h <- est.h(rt_ru.ror2.dat, loci.data = rt.ror2.loci.dat, fixed = FALSE)
rt_ty.ror2.h <- est.h(rt_ty.ror2.dat, loci.data = rt.ror2.loci.dat, fixed = FALSE)
rt_rt.ror2.h <- est.h(rt_rt.ror2.dat, loci.data = rt.ror2.loci.dat, fixed = FALSE)

rt_ru.apc.camk4.h <- est.h(rt_ru.apc.camk4.dat, loci.data = rt.apc.camk4.loci.dat, fixed = FALSE)
rt_ty.apc.camk4.h <- est.h(rt_ty.apc.camk4.dat, loci.data = rt.apc.camk4.loci.dat, fixed = FALSE)
rt_rt.apc.camk4.h <- est.h(rt_rt.apc.camk4.dat, loci.data = rt.apc.camk4.loci.dat, fixed = FALSE)

rt_ru.ice1.lncrna.h <- est.h(rt_ru.ice1.lncrna.dat, loci.data = rt.ice1.lncrna.loci.dat, fixed = FALSE)
rt_ty.ice1.lncrna.h <- est.h(rt_ty.ice1.lncrna.dat, loci.data = rt.ice1.lncrna.loci.dat, fixed = FALSE)
rt_rt.ice1.lncrna.h <- est.h(rt_rt.ice1.lncrna.dat, loci.data = rt.ice1.lncrna.loci.dat, fixed = FALSE)

rt_ru.pde1c.h <- est.h(rt_ru.pde1c.dat, loci.data = rt.pde1c.loci.dat, fixed = FALSE)
rt_ty.pde1c.h <- est.h(rt_ty.pde1c.dat, loci.data = rt.pde1c.loci.dat, fixed = FALSE)
rt_rt.pde1c.h <- est.h(rt_rt.pde1c.dat, loci.data = rt.pde1c.loci.dat, fixed = FALSE)

## Query hybrid index values and write to `data.csv`
print(rt_ru.kitlg.h)
print(rt_rt.kitlg.h)
print(rt_ty.kitlg.h)

print(rt_ru.plxnc1.h)
print(rt_rt.plxnc1.h)
print(rt_ty.plxnc1.h)

print(rt_ru.spef2.prlr.h)
print(rt_rt.spef2.prlr.h)
print(rt_ty.spef2.prlr.h)

print(rt_ru.slc45a2.h)
print(rt_rt.slc45a2.h)
print(rt_ty.slc45a2.h)

print(rt_ru.tyrp1.h)
print(rt_rt.tyrp1.h)
print(rt_ty.tyrp1.h)

print(rt_ru.bnc2.h)
print(rt_rt.bnc2.h)
print(rt_ty.bnc2.h)

print(rt_ru.bnc2.gwa.h)
print(rt_rt.bnc2.gwa.h)
print(rt_ty.bnc2.gwa.h)

print(rt_ru.gnaq.h)
print(rt_rt.gnaq.h)
print(rt_ty.gnaq.h)

print(rt_ru.ror2.h)
print(rt_rt.ror2.h)
print(rt_ty.ror2.h)

print(rt_ru.apc.camk4.h)
print(rt_rt.apc.camk4.h)
print(rt_ty.apc.camk4.h)

print(rt_ru.ice1.lncrna.h)
print(rt_rt.ice1.lncrna.h)
print(rt_ty.ice1.lncrna.h)

print(rt_ru.pde1c.h)
print(rt_rt.pde1c.h)
print(rt_ty.pde1c.h)

### Read in and format data for HZAR (rustica-tytleri)----------------------

# These data include the sample, population, lat, long, locality, hybrid index
# (background), and hybrid index (outlier regions).

data.rt <- read.csv('./data.rustica-tytleri.csv',header=T)
names(data.rt)

## Standardize phenotype vector ranges between 0 and 1 and add to data frame
mean_rwl_norm = (data.rt$mean_rwl-min(data.rt$mean_rwl,na.rm = T))/(max(data.rt$mean_rwl,na.rm = T)-min(data.rt$mean_rwl,na.rm = T))
mean_ts_random_norm = (data.rt$mean_ts_random-min(data.rt$mean_ts_random,na.rm = T))/(max(data.rt$mean_ts_random,na.rm = T)-min(data.rt$mean_ts_random,na.rm = T))
breast_avg_bright_norm = (data.rt$breast_avg_bright-min(data.rt$breast_avg_bright,na.rm = T))/(max(data.rt$breast_avg_bright,na.rm = T)-min(data.rt$breast_avg_bright,na.rm = T))

data.rt$wing <- mean_rwl_norm
data.rt$tail <- mean_ts_random_norm
data.rt$vent <- breast_avg_bright_norm

## Subset to get locality, long, lat, HI-background, and HI-outlier

data.rt <- dplyr::select(data.rt, loc,long, lat, hi_kitlg, hi_plxnc1, hi_spef2_prlr, hi_slc45a2, hi_tyrp1, hi_bnc2, hi_gnaq, hi_ror2, hi_apc_camk4, hi_ice1_lncrna, hi_pde1c)
head(data.rt)

## Calculate means, variances, and counts
data.rt.melt<-reshape2::melt(data.rt, id.vars=c("loc"))
str(data.rt.melt)
head(data.rt.melt)
data.rt.hzar<- data.rt.melt %>%
  group_by(loc, variable) %>%
  dplyr::summarise(mean=mean(value, na.rm=TRUE), var=var(value, na.rm=TRUE), count=sum(!is.na(value)))
head(data.rt.hzar)

## Change to wide format
d<-reshape2::melt(data.rt.hzar, id.vars=c("loc", "variable"))
head(d)
colnames(d)[3]<-"metric"
data.rt.hzar.cast<-dcast(d, loc~ variable + metric )
print(data.rt.hzar.cast)

## Calculate distances (km) for transect
rt.dists<-hzar.map.latLongSites(data.rt.hzar.cast$loc, data.rt.hzar.cast$lat_mean, data.rt.hzar.cast$long_mean, degrees = TRUE)
join_func<-function(y){
  y<-left_join(y, data.rt.hzar.cast, by=c("site"="loc"))
  y<-y[,-c(2:3)]
  y<-arrange(y, km)
}

## Subset points to the sampling transects
rt.dists<-filter(rt.dists, long.deg>70)
rt.dists<-filter(rt.dists, !site=="urumqi") #urumqi the straight line
rt.dists<-filter(rt.dists, !site=="krasnoyarsk2") #urumqi the straight line
print(rt.dists)

##Add km transect distances
#rt.dists$km <- hzar.map.distanceFromSite(rt.dists,"yekaterinburg",units="Km")
rt.dists$km <- hzar.map.distanceFromSite(rt.dists,"karasuk",units="Km")
#rt.dists$km <- hzar.map.distanceFromSite(rt.dists,"krasnoyarsk",units="Km")
rt.dists<-join_func(rt.dists)
rt.dists<-na.omit(rt.dists)
print(rt.dists)

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
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$kitlg <- list();
## Space to hold the observed data
rt$kitlg$obs <- list();
## Space to hold the models to fit
rt$kitlg$models <- list();
## Space to hold the compiled fit requests
rt$kitlg$fitRs <- list();
## Space to hold the output data chains
rt$kitlg$runs <- list();
## Space to hold the analysed data
rt$kitlg$analysis <- list();

## Assign data for KITLG region
rt$kitlg$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_kitlg_mean,
                             rt.dists$hi_kitlg_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$kitlg$obs);

## Make a helper function to set cline models
rt.loadkitlgmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  rt$kitlg$models[[id]] <<- hzar.makeCline1DFreq(rt$kitlg$obs, scaling, tails)

rt.loadkitlgmodel("fixed","none","modelI");
rt.loadkitlgmodel("free" ,"none","modelII");
rt.loadkitlgmodel("free" ,"both","modelIII");
rt.loadkitlgmodel("free" ,"right","modelIV");
rt.loadkitlgmodel("free" ,"left","modelV");
rt.loadkitlgmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rt$kitlg$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$kitlg$models <- sapply(rt$kitlg$models,
                          hzar.model.addBoxReq,
                          -30 , 1973,
                          simplify=FALSE)

#rt$kitlg$models <- sapply(rt$kitlg$models,
#                          hzar.model.addBoxReq,
#                          -30 , 1383,
#                          simplify=FALSE)

## Check the updated settings
print(rt$kitlg$models)

## Compile each of the models to prepare for fitting
rt$kitlg$fitRs$init <- sapply(rt$kitlg$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=rt$kitlg$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
rt$kitlg$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$kitlg$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$kitlg$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$kitlg$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$kitlg$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$kitlg$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$kitlg$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$kitlg$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$kitlg$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$kitlg$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$kitlg$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$kitlg$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$kitlg$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$kitlg$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$kitlg$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$kitlg$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$kitlg$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$kitlg$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$kitlg$fitRs$init)

## Run the models for an initial chain
rt$kitlg$runs$init <- list()

rt$kitlg$runs$init$modelI <-
  hzar.doFit(rt$kitlg$fitRs$init$modelI)

rt$kitlg$runs$init$modelII <-
  hzar.doFit(rt$kitlg$fitRs$init$modelII)

rt$kitlg$runs$init$modelIII <-
  hzar.doFit(rt$kitlg$fitRs$init$modelIII)

rt$kitlg$runs$init$modelIV <-
  hzar.doFit(rt$kitlg$fitRs$init$modelIV)

rt$kitlg$runs$init$modelV <-
  hzar.doFit(rt$kitlg$fitRs$init$modelV)

rt$kitlg$runs$init$modelVI <-
  hzar.doFit(rt$kitlg$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$kitlg$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$kitlg$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$kitlg$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$kitlg$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$kitlg$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$kitlg$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$kitlg$fitRs$chains <-
  lapply(rt$kitlg$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$kitlg$fitRs$chains <-
  hzar.multiFitRequest(rt$kitlg$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$kitlg$runs$chains <-  hzar.doChain.multi(rt$kitlg$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$kitlg$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$kitlg$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$kitlg$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$kitlg$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$kitlg$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$kitlg$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$kitlg$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$kitlg$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$kitlg$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$kitlg$runs$init$modelI)
rt$kitlg$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$kitlg$runs$init$modelII)
rt$kitlg$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$kitlg$runs$init$modelIII)
rt$kitlg$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$kitlg$runs$init$modelIV)
rt$kitlg$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$kitlg$runs$init$modelV)
rt$kitlg$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$kitlg$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$kitlg$analysis$oDG <-
  hzar.make.obsDataGroup(rt$kitlg$analysis$initDGs)
rt$kitlg$analysis$oDG <-
  hzar.copyModelLabels(rt$kitlg$analysis$initDGs,
                       rt$kitlg$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$kitlg$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$kitlg$runs$chains,
                                hzar.dataGroup.add),
                         rt$kitlg$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$kitlg$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$kitlg$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$kitlg$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$kitlg$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$kitlg$analysis$model.name <-
        rownames(rt$kitlg$analysis$AICcTable
        )[[ which.min(rt$kitlg$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$kitlg$analysis$model.selected <-
  rt$kitlg$analysis$oDG$data.groups[[rt$kitlg$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$kitlg$analysis$model.selected,
                         names(rt$kitlg$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$kitlg$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$kitlg$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$kitlg$analysis$model.selected);

# End Analysis

#dev.off()

### Analysis on PLXNC1 region (rustica-tytleri)-----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$plxnc1 <- list();
## Space to hold the observed data
rt$plxnc1$obs <- list();
## Space to hold the models to fit
rt$plxnc1$models <- list();
## Space to hold the compiled fit requests
rt$plxnc1$fitRs <- list();
## Space to hold the output data chains
rt$plxnc1$runs <- list();
## Space to hold the analysed data
rt$plxnc1$analysis <- list();

## Assign data for PLXNC1 region
rt$plxnc1$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_plxnc1_mean,
                             rt.dists$hi_plxnc1_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$plxnc1$obs);

## Make a helper function to set cline models
rt.loadplxnc1model <- function(scaling,tails,
                               id=paste(scaling,tails,sep="."))
  rt$plxnc1$models[[id]] <<- hzar.makeCline1DFreq(rt$plxnc1$obs, scaling, tails)

rt.loadplxnc1model("fixed","none","modelI");
rt.loadplxnc1model("free" ,"none","modelII");
rt.loadplxnc1model("free" ,"both","modelIII");
rt.loadplxnc1model("free" ,"right","modelIV");
rt.loadplxnc1model("free" ,"left","modelV");
rt.loadplxnc1model("free" ,"mirror","modelVI");

## Check the default settings
print(rt$plxnc1$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$plxnc1$models <- sapply(rt$plxnc1$models,
                           hzar.model.addBoxReq,
                           -30 , 1973,
                           simplify=FALSE)

## Check the updated settings
print(rt$plxnc1$models)

## Compile each of the models to prepare for fitting
rt$plxnc1$fitRs$init <- sapply(rt$plxnc1$models,
                               hzar.first.fitRequest.old.ML,
                               obsData=rt$plxnc1$obs,
                               verbose=FALSE,
                               simplify=FALSE)

## Update the settings for the fitter if desired.
rt$plxnc1$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$plxnc1$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$plxnc1$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$plxnc1$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$plxnc1$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$plxnc1$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$plxnc1$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$plxnc1$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$plxnc1$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$plxnc1$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$plxnc1$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$plxnc1$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$plxnc1$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$plxnc1$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$plxnc1$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$plxnc1$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$plxnc1$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$plxnc1$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$plxnc1$fitRs$init)

## Run the models for an initial chain
rt$plxnc1$runs$init <- list()

rt$plxnc1$runs$init$modelI <-
  hzar.doFit(rt$plxnc1$fitRs$init$modelI)

rt$plxnc1$runs$init$modelII <-
  hzar.doFit(rt$plxnc1$fitRs$init$modelII)

rt$plxnc1$runs$init$modelIII <-
  hzar.doFit(rt$plxnc1$fitRs$init$modelIII)

rt$plxnc1$runs$init$modelIV <-
  hzar.doFit(rt$plxnc1$fitRs$init$modelIV)

rt$plxnc1$runs$init$modelV <-
  hzar.doFit(rt$plxnc1$fitRs$init$modelV)

rt$plxnc1$runs$init$modelVI <-
  hzar.doFit(rt$plxnc1$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$plxnc1$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$plxnc1$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$plxnc1$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$plxnc1$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$plxnc1$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$plxnc1$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$plxnc1$fitRs$chains <-
  lapply(rt$plxnc1$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$plxnc1$fitRs$chains <-
  hzar.multiFitRequest(rt$plxnc1$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$plxnc1$runs$chains <-  hzar.doChain.multi(rt$plxnc1$fitRs$chains,
                                             doPar=TRUE,
                                             inOrder=FALSE,
                                             count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$plxnc1$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$plxnc1$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$plxnc1$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$plxnc1$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$plxnc1$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$plxnc1$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$plxnc1$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$plxnc1$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$plxnc1$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$plxnc1$runs$init$modelI)
rt$plxnc1$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$plxnc1$runs$init$modelII)
rt$plxnc1$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$plxnc1$runs$init$modelIII)
rt$plxnc1$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$plxnc1$runs$init$modelIV)
rt$plxnc1$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$plxnc1$runs$init$modelV)
rt$plxnc1$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$plxnc1$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$plxnc1$analysis$oDG <-
  hzar.make.obsDataGroup(rt$plxnc1$analysis$initDGs)
rt$plxnc1$analysis$oDG <-
  hzar.copyModelLabels(rt$plxnc1$analysis$initDGs,
                       rt$plxnc1$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$plxnc1$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$plxnc1$runs$chains,
                                hzar.dataGroup.add),
                         rt$plxnc1$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$plxnc1$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$plxnc1$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$plxnc1$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$plxnc1$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$plxnc1$analysis$model.name <-
        rownames(rt$plxnc1$analysis$AICcTable
        )[[ which.min(rt$plxnc1$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$plxnc1$analysis$model.selected <-
  rt$plxnc1$analysis$oDG$data.groups[[rt$plxnc1$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$plxnc1$analysis$model.selected,
                         names(rt$plxnc1$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$plxnc1$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$plxnc1$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$plxnc1$analysis$model.selected);
rect(xleft=rt$plxnc1$analysis$model.selected$ML.cline$param.all$center-((rt$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=rt$plxnc1$analysis$model.selected$ML.cline$param.all$center+((rt$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#E28026')

## End Analysis

#dev.off()

### Analysis on SPEF2/PRLR region (rustica-tytleri)-------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$spef2.prlr <- list();
## Space to hold the observed data
rt$spef2.prlr$obs <- list();
## Space to hold the models to fit
rt$spef2.prlr$models <- list();
## Space to hold the compiled fit requests
rt$spef2.prlr$fitRs <- list();
## Space to hold the output data chains
rt$spef2.prlr$runs <- list();
## Space to hold the analysed data
rt$spef2.prlr$analysis <- list();

## Assign data for SPEF2/PRLR region
rt$spef2.prlr$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_spef2_prlr_mean,
                             rt.dists$hi_spef2_prlr_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$spef2.prlr$obs);

## Make a helper function to set cline models
rt.loadspef2.prlrmodel <- function(scaling,tails,
                                   id=paste(scaling,tails,sep="."))
  rt$spef2.prlr$models[[id]] <<- hzar.makeCline1DFreq(rt$spef2.prlr$obs, scaling, tails)

rt.loadspef2.prlrmodel("fixed","none","modelI");
rt.loadspef2.prlrmodel("free" ,"none","modelII");
rt.loadspef2.prlrmodel("free" ,"both","modelIII");
rt.loadspef2.prlrmodel("free" ,"right","modelIV");
rt.loadspef2.prlrmodel("free" ,"left","modelV");
rt.loadspef2.prlrmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rt$spef2.prlr$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$spef2.prlr$models <- sapply(rt$spef2.prlr$models,
                               hzar.model.addBoxReq,
                               -30 , 1973,
                               simplify=FALSE)

## Check the updated settings
print(rt$spef2.prlr$models)

## Compile each of the models to prepare for fitting
rt$spef2.prlr$fitRs$init <- sapply(rt$spef2.prlr$models,
                                   hzar.first.fitRequest.old.ML,
                                   obsData=rt$spef2.prlr$obs,
                                   verbose=FALSE,
                                   simplify=FALSE)

## Update the settings for the fitter if desired.
rt$spef2.prlr$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$spef2.prlr$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$spef2.prlr$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$spef2.prlr$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$spef2.prlr$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$spef2.prlr$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$spef2.prlr$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$spef2.prlr$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$spef2.prlr$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$spef2.prlr$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$spef2.prlr$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$spef2.prlr$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$spef2.prlr$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$spef2.prlr$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$spef2.prlr$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$spef2.prlr$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$spef2.prlr$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$spef2.prlr$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$spef2.prlr$fitRs$init)

## Run the models for an initial chain
rt$spef2.prlr$runs$init <- list()

rt$spef2.prlr$runs$init$modelI <-
  hzar.doFit(rt$spef2.prlr$fitRs$init$modelI)

rt$spef2.prlr$runs$init$modelII <-
  hzar.doFit(rt$spef2.prlr$fitRs$init$modelII)

rt$spef2.prlr$runs$init$modelIII <-
  hzar.doFit(rt$spef2.prlr$fitRs$init$modelIII)

rt$spef2.prlr$runs$init$modelIV <-
  hzar.doFit(rt$spef2.prlr$fitRs$init$modelIV)

rt$spef2.prlr$runs$init$modelV <-
  hzar.doFit(rt$spef2.prlr$fitRs$init$modelV)

rt$spef2.prlr$runs$init$modelVI <-
  hzar.doFit(rt$spef2.prlr$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$spef2.prlr$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$spef2.prlr$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$spef2.prlr$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$spef2.prlr$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$spef2.prlr$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$spef2.prlr$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$spef2.prlr$fitRs$chains <-
  lapply(rt$spef2.prlr$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$spef2.prlr$fitRs$chains <-
  hzar.multiFitRequest(rt$spef2.prlr$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$spef2.prlr$runs$chains <-  hzar.doChain.multi(rt$spef2.prlr$fitRs$chains,
                                                 doPar=TRUE,
                                                 inOrder=FALSE,
                                                 count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$spef2.prlr$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$spef2.prlr$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$spef2.prlr$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$spef2.prlr$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$spef2.prlr$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$spef2.prlr$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$spef2.prlr$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$spef2.prlr$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$spef2.prlr$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$spef2.prlr$runs$init$modelI)
rt$spef2.prlr$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$spef2.prlr$runs$init$modelII)
rt$spef2.prlr$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$spef2.prlr$runs$init$modelIII)
rt$spef2.prlr$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$spef2.prlr$runs$init$modelIV)
rt$spef2.prlr$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$spef2.prlr$runs$init$modelV)
rt$spef2.prlr$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$spef2.prlr$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$spef2.prlr$analysis$oDG <-
  hzar.make.obsDataGroup(rt$spef2.prlr$analysis$initDGs)
rt$spef2.prlr$analysis$oDG <-
  hzar.copyModelLabels(rt$spef2.prlr$analysis$initDGs,
                       rt$spef2.prlr$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$spef2.prlr$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$spef2.prlr$runs$chains,
                                hzar.dataGroup.add),
                         rt$spef2.prlr$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$spef2.prlr$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$spef2.prlr$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$spef2.prlr$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$spef2.prlr$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$spef2.prlr$analysis$model.name <-
        rownames(rt$spef2.prlr$analysis$AICcTable
        )[[ which.min(rt$spef2.prlr$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$spef2.prlr$analysis$model.selected <-
  rt$spef2.prlr$analysis$oDG$data.groups[[rt$spef2.prlr$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$spef2.prlr$analysis$model.selected,
                         names(rt$spef2.prlr$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$spef2.prlr$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$spef2.prlr$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$spef2.prlr$analysis$model.selected);
rect(xleft=rt$spef2.prlr$analysis$model.selected$ML.cline$param.all$center-((rt$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)/2),ytop=0,xright=rt$spef2.prlr$analysis$model.selected$ML.cline$param.all$center+((rt$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)/2),ybottom=-0.05,col='#E28026')

## Plot the 95% credible cline for model II (free; none) --> very similar support to model I (delta AIC = 0.04)
hzar.plot.fzCline(rt$spef2.prlr$analysis$oDG$data.groups$modelII);
rect(xleft=rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),ytop=0,xright=rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),ybottom=-0.05,col='#E28026')

## End Analysis

#dev.off()

### Analysis on SLC45A2 region (rustica-tytleri)----------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$slc45a2 <- list();
## Space to hold the observed data
rt$slc45a2$obs <- list();
## Space to hold the models to fit
rt$slc45a2$models <- list();
## Space to hold the compiled fit requests
rt$slc45a2$fitRs <- list();
## Space to hold the output data chains
rt$slc45a2$runs <- list();
## Space to hold the analysed data
rt$slc45a2$analysis <- list();

## Assign data for SLC45A2 region
rt$slc45a2$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_slc45a2_mean,
                             rt.dists$hi_slc45a2_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$slc45a2$obs);

## Make a helper function to set cline models
rt.loadslc45a2model <- function(scaling,tails,
                                id=paste(scaling,tails,sep="."))
  rt$slc45a2$models[[id]] <<- hzar.makeCline1DFreq(rt$slc45a2$obs, scaling, tails)

rt.loadslc45a2model("fixed","none","modelI");
rt.loadslc45a2model("free" ,"none","modelII");
rt.loadslc45a2model("free" ,"both","modelIII");
rt.loadslc45a2model("free" ,"right","modelIV");
rt.loadslc45a2model("free" ,"left","modelV");
rt.loadslc45a2model("free" ,"mirror","modelVI");

## Check the default settings
print(rt$slc45a2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$slc45a2$models <- sapply(rt$slc45a2$models,
                            hzar.model.addBoxReq,
                            -30 , 1973,
                            simplify=FALSE)

## Check the updated settings
print(rt$slc45a2$models)

## Compile each of the models to prepare for fitting
rt$slc45a2$fitRs$init <- sapply(rt$slc45a2$models,
                                hzar.first.fitRequest.old.ML,
                                obsData=rt$slc45a2$obs,
                                verbose=FALSE,
                                simplify=FALSE)

## Update the settings for the fitter if desired.
rt$slc45a2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$slc45a2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$slc45a2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$slc45a2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$slc45a2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$slc45a2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$slc45a2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$slc45a2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$slc45a2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$slc45a2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$slc45a2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$slc45a2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$slc45a2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$slc45a2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$slc45a2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$slc45a2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$slc45a2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$slc45a2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$slc45a2$fitRs$init)

## Run the models for an initial chain
rt$slc45a2$runs$init <- list()

rt$slc45a2$runs$init$modelI <-
  hzar.doFit(rt$slc45a2$fitRs$init$modelI)

rt$slc45a2$runs$init$modelII <-
  hzar.doFit(rt$slc45a2$fitRs$init$modelII)

rt$slc45a2$runs$init$modelIII <-
  hzar.doFit(rt$slc45a2$fitRs$init$modelIII)

rt$slc45a2$runs$init$modelIV <-
  hzar.doFit(rt$slc45a2$fitRs$init$modelIV)

rt$slc45a2$runs$init$modelV <-
  hzar.doFit(rt$slc45a2$fitRs$init$modelV)

rt$slc45a2$runs$init$modelVI <-
  hzar.doFit(rt$slc45a2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$slc45a2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$slc45a2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$slc45a2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$slc45a2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$slc45a2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$slc45a2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$slc45a2$fitRs$chains <-
  lapply(rt$slc45a2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$slc45a2$fitRs$chains <-
  hzar.multiFitRequest(rt$slc45a2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$slc45a2$runs$chains <-  hzar.doChain.multi(rt$slc45a2$fitRs$chains,
                                              doPar=TRUE,
                                              inOrder=FALSE,
                                              count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$slc45a2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$slc45a2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$slc45a2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$slc45a2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$slc45a2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$slc45a2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$slc45a2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$slc45a2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$slc45a2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$slc45a2$runs$init$modelI)
rt$slc45a2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$slc45a2$runs$init$modelII)
rt$slc45a2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$slc45a2$runs$init$modelIII)
rt$slc45a2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$slc45a2$runs$init$modelIV)
rt$slc45a2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$slc45a2$runs$init$modelV)
rt$slc45a2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$slc45a2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$slc45a2$analysis$oDG <-
  hzar.make.obsDataGroup(rt$slc45a2$analysis$initDGs)
rt$slc45a2$analysis$oDG <-
  hzar.copyModelLabels(rt$slc45a2$analysis$initDGs,
                       rt$slc45a2$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$slc45a2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$slc45a2$runs$chains,
                                hzar.dataGroup.add),
                         rt$slc45a2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$slc45a2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$slc45a2$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$slc45a2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$slc45a2$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$slc45a2$analysis$model.name <-
        rownames(rt$slc45a2$analysis$AICcTable
        )[[ which.min(rt$slc45a2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$slc45a2$analysis$model.selected <-
  rt$slc45a2$analysis$oDG$data.groups[[rt$slc45a2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$slc45a2$analysis$model.selected,
                         names(rt$slc45a2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$slc45a2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$slc45a2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$slc45a2$analysis$model.selected);

## Plot the 95% credible cline for model II (free; none) --> very similar support to model I (delta AIC = 1.73)
hzar.plot.fzCline(rt$slc45a2$analysis$oDG$data.groups$modelII);
rect(xleft=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),ytop=0,xright=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),ybottom=-0.05,col='#E28026')

## End Analysis

#dev.off()

### Analysis on BNC2 region (rustica-tytleri)-------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$bnc2 <- list();
## Space to hold the observed data
rt$bnc2$obs <- list();
## Space to hold the models to fit
rt$bnc2$models <- list();
## Space to hold the compiled fit requests
rt$bnc2$fitRs <- list();
## Space to hold the output data chains
rt$bnc2$runs <- list();
## Space to hold the analysed data
rt$bnc2$analysis <- list();

## Assign data for BNC2 region
rt$bnc2$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_bnc2_mean,
                             rt.dists$hi_bnc2_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$bnc2$obs);

## Make a helper function to set cline models
rt.loadbnc2model <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  rt$bnc2$models[[id]] <<- hzar.makeCline1DFreq(rt$bnc2$obs, scaling, tails)

rt.loadbnc2model("fixed","none","modelI");
rt.loadbnc2model("free" ,"none","modelII");
rt.loadbnc2model("free" ,"both","modelIII");
rt.loadbnc2model("free" ,"right","modelIV");
rt.loadbnc2model("free" ,"left","modelV");
rt.loadbnc2model("free" ,"mirror","modelVI");

## Check the default settings
print(rt$bnc2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$bnc2$models <- sapply(rt$bnc2$models,
                         hzar.model.addBoxReq,
                         -30 , 1973,
                         simplify=FALSE)

## Check the updated settings
print(rt$bnc2$models)

## Compile each of the models to prepare for fitting
rt$bnc2$fitRs$init <- sapply(rt$bnc2$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=rt$bnc2$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
rt$bnc2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$bnc2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$bnc2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$bnc2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$bnc2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$bnc2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$bnc2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$bnc2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$bnc2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$bnc2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$bnc2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$bnc2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$bnc2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$bnc2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$bnc2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$bnc2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$bnc2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$bnc2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$bnc2$fitRs$init)

## Run the models for an initial chain
rt$bnc2$runs$init <- list()

rt$bnc2$runs$init$modelI <-
  hzar.doFit(rt$bnc2$fitRs$init$modelI)

rt$bnc2$runs$init$modelII <-
  hzar.doFit(rt$bnc2$fitRs$init$modelII)

rt$bnc2$runs$init$modelIII <-
  hzar.doFit(rt$bnc2$fitRs$init$modelIII)

rt$bnc2$runs$init$modelIV <-
  hzar.doFit(rt$bnc2$fitRs$init$modelIV)

rt$bnc2$runs$init$modelV <-
  hzar.doFit(rt$bnc2$fitRs$init$modelV)

rt$bnc2$runs$init$modelVI <-
  hzar.doFit(rt$bnc2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$bnc2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$bnc2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$bnc2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$bnc2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$bnc2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$bnc2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$bnc2$fitRs$chains <-
  lapply(rt$bnc2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$bnc2$fitRs$chains <-
  hzar.multiFitRequest(rt$bnc2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$bnc2$runs$chains <-  hzar.doChain.multi(rt$bnc2$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$bnc2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$bnc2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$bnc2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$bnc2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$bnc2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$bnc2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$bnc2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$bnc2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$bnc2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$bnc2$runs$init$modelI)
rt$bnc2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$bnc2$runs$init$modelII)
rt$bnc2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$bnc2$runs$init$modelIII)
rt$bnc2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$bnc2$runs$init$modelIV)
rt$bnc2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$bnc2$runs$init$modelV)
rt$bnc2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$bnc2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$bnc2$analysis$oDG <-
  hzar.make.obsDataGroup(rt$bnc2$analysis$initDGs)
rt$bnc2$analysis$oDG <-
  hzar.copyModelLabels(rt$bnc2$analysis$initDGs,
                       rt$bnc2$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$bnc2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$bnc2$runs$chains,
                                hzar.dataGroup.add),
                         rt$bnc2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$bnc2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$bnc2$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$bnc2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$bnc2$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$bnc2$analysis$model.name <-
        rownames(rt$bnc2$analysis$AICcTable
        )[[ which.min(rt$bnc2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$bnc2$analysis$model.selected <-
  rt$bnc2$analysis$oDG$data.groups[[rt$bnc2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$bnc2$analysis$model.selected,
                         names(rt$bnc2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$bnc2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$bnc2$analysis$model.selected);
hzar.plot.cline(rt$bnc2$analysis$initDGs$modelII);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$bnc2$analysis$model.selected);
rect(xleft=rt$bnc2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt$bnc2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),ytop=0,xright=rt$bnc2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt$bnc2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),ybottom=-0.05,col='#E28026')
rect(xleft=rt$bnc2$analysis$oDG$data.groups$modelI$ML.cline$param.all$center-((rt$bnc2$analysis$oDG$data.groups$modelI$ML.cline$param.all$width)/2),ytop=0,xright=rt$bnc2$analysis$oDG$data.groups$modelI$ML.cline$param.all$center+((rt$bnc2$analysis$oDG$data.groups$modelI$ML.cline$param.all$width)/2),ybottom=-0.05,col='#E28026')
hzar.plot.fzCline(rt$bnc2$analysis$initDGs$modelII);

## End Analysis

#dev.off()

### Analysis on GNAQ region (rustica-tytleri)-------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$gnaq <- list();
## Space to hold the observed data
rt$gnaq$obs <- list();
## Space to hold the models to fit
rt$gnaq$models <- list();
## Space to hold the compiled fit requests
rt$gnaq$fitRs <- list();
## Space to hold the output data chains
rt$gnaq$runs <- list();
## Space to hold the analysed data
rt$gnaq$analysis <- list();

## Assign data for GNAQ region
rt$gnaq$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_gnaq_mean,
                             rt.dists$hi_gnaq_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$gnaq$obs);

## Make a helper function to set cline models
rt.loadgnaqmodel <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  rt$gnaq$models[[id]] <<- hzar.makeCline1DFreq(rt$gnaq$obs, scaling, tails)

rt.loadgnaqmodel("fixed","none","modelI");
rt.loadgnaqmodel("free" ,"none","modelII");
rt.loadgnaqmodel("free" ,"both","modelIII");
rt.loadgnaqmodel("free" ,"right","modelIV");
rt.loadgnaqmodel("free" ,"left","modelV");
rt.loadgnaqmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rt$gnaq$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$gnaq$models <- sapply(rt$gnaq$models,
                         hzar.model.addBoxReq,
                         -30 , 1973,
                         simplify=FALSE)

## Check the updated settings
print(rt$gnaq$models)

## Compile each of the models to prepare for fitting
rt$gnaq$fitRs$init <- sapply(rt$gnaq$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=rt$gnaq$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
rt$gnaq$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$gnaq$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$gnaq$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$gnaq$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$gnaq$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$gnaq$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$gnaq$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$gnaq$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$gnaq$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$gnaq$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$gnaq$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$gnaq$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$gnaq$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$gnaq$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$gnaq$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$gnaq$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$gnaq$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$gnaq$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$gnaq$fitRs$init)

## Run the models for an initial chain
rt$gnaq$runs$init <- list()

rt$gnaq$runs$init$modelI <-
  hzar.doFit(rt$gnaq$fitRs$init$modelI)

rt$gnaq$runs$init$modelII <-
  hzar.doFit(rt$gnaq$fitRs$init$modelII)

rt$gnaq$runs$init$modelIII <-
  hzar.doFit(rt$gnaq$fitRs$init$modelIII)

rt$gnaq$runs$init$modelIV <-
  hzar.doFit(rt$gnaq$fitRs$init$modelIV)

rt$gnaq$runs$init$modelV <-
  hzar.doFit(rt$gnaq$fitRs$init$modelV)

rt$gnaq$runs$init$modelVI <-
  hzar.doFit(rt$gnaq$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$gnaq$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$gnaq$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$gnaq$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$gnaq$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$gnaq$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$gnaq$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$gnaq$fitRs$chains <-
  lapply(rt$gnaq$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$gnaq$fitRs$chains <-
  hzar.multiFitRequest(rt$gnaq$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$gnaq$runs$chains <-  hzar.doChain.multi(rt$gnaq$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$gnaq$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$gnaq$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$gnaq$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$gnaq$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$gnaq$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$gnaq$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$gnaq$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$gnaq$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$gnaq$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$gnaq$runs$init$modelI)
rt$gnaq$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$gnaq$runs$init$modelII)
rt$gnaq$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$gnaq$runs$init$modelIII)
rt$gnaq$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$gnaq$runs$init$modelIV)
rt$gnaq$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$gnaq$runs$init$modelV)
rt$gnaq$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$gnaq$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$gnaq$analysis$oDG <-
  hzar.make.obsDataGroup(rt$gnaq$analysis$initDGs)
rt$gnaq$analysis$oDG <-
  hzar.copyModelLabels(rt$gnaq$analysis$initDGs,
                       rt$gnaq$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$gnaq$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$gnaq$runs$chains,
                                hzar.dataGroup.add),
                         rt$gnaq$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$gnaq$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$gnaq$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$gnaq$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$gnaq$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$gnaq$analysis$model.name <-
        rownames(rt$gnaq$analysis$AICcTable
        )[[ which.min(rt$gnaq$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$gnaq$analysis$model.selected <-
  rt$gnaq$analysis$oDG$data.groups[[rt$gnaq$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$gnaq$analysis$model.selected,
                         names(rt$gnaq$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$gnaq$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$gnaq$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$gnaq$analysis$model.selected);

## End Analysis

#dev.off()

par(mfrow=c(4,1))
hzar.plot.fzCline(rt$back$analysis$model.selected);
hzar.plot.fzCline(rt$kitlg$analysis$model.selected);
hzar.plot.fzCline(rt$bnc2$analysis$model.selected);
hzar.plot.fzCline(rt$apc.camk4$analysis$model.selected);

### Analysis on ROR2 region (rustica-tytleri)-------------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$ror2 <- list();
## Space to hold the observed data
rt$ror2$obs <- list();
## Space to hold the models to fit
rt$ror2$models <- list();
## Space to hold the compiled fit requests
rt$ror2$fitRs <- list();
## Space to hold the output data chains
rt$ror2$runs <- list();
## Space to hold the analysed data
rt$ror2$analysis <- list();

## Assign data for ROR2 region
rt$ror2$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_ror2_mean,
                             rt.dists$hi_ror2_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$ror2$obs);

## Make a helper function to set cline models
rt.loadror2model <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  rt$ror2$models[[id]] <<- hzar.makeCline1DFreq(rt$ror2$obs, scaling, tails)

rt.loadror2model("fixed","none","modelI");
rt.loadror2model("free" ,"none","modelII");
rt.loadror2model("free" ,"both","modelIII");
rt.loadror2model("free" ,"right","modelIV");
rt.loadror2model("free" ,"left","modelV");
rt.loadror2model("free" ,"mirror","modelVI");

## Check the default settings
print(rt$ror2$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$ror2$models <- sapply(rt$ror2$models,
                         hzar.model.addBoxReq,
                         -30 , 1973,
                         simplify=FALSE)

## Check the updated settings
print(rt$ror2$models)

## Compile each of the models to prepare for fitting
rt$ror2$fitRs$init <- sapply(rt$ror2$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=rt$ror2$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
rt$ror2$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ror2$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ror2$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$ror2$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ror2$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ror2$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$ror2$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ror2$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ror2$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$ror2$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ror2$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ror2$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$ror2$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ror2$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ror2$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$ror2$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ror2$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ror2$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$ror2$fitRs$init)

## Run the models for an initial chain
rt$ror2$runs$init <- list()

rt$ror2$runs$init$modelI <-
  hzar.doFit(rt$ror2$fitRs$init$modelI)

rt$ror2$runs$init$modelII <-
  hzar.doFit(rt$ror2$fitRs$init$modelII)

rt$ror2$runs$init$modelIII <-
  hzar.doFit(rt$ror2$fitRs$init$modelIII)

rt$ror2$runs$init$modelIV <-
  hzar.doFit(rt$ror2$fitRs$init$modelIV)

rt$ror2$runs$init$modelV <-
  hzar.doFit(rt$ror2$fitRs$init$modelV)

rt$ror2$runs$init$modelVI <-
  hzar.doFit(rt$ror2$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$ror2$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$ror2$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$ror2$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$ror2$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$ror2$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$ror2$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$ror2$fitRs$chains <-
  lapply(rt$ror2$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$ror2$fitRs$chains <-
  hzar.multiFitRequest(rt$ror2$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$ror2$runs$chains <-  hzar.doChain.multi(rt$ror2$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$ror2$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$ror2$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$ror2$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$ror2$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$ror2$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$ror2$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$ror2$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$ror2$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$ror2$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$ror2$runs$init$modelI)
rt$ror2$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$ror2$runs$init$modelII)
rt$ror2$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$ror2$runs$init$modelIII)
rt$ror2$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$ror2$runs$init$modelIV)
rt$ror2$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$ror2$runs$init$modelV)
rt$ror2$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$ror2$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$ror2$analysis$oDG <-
  hzar.make.obsDataGroup(rt$ror2$analysis$initDGs)
rt$ror2$analysis$oDG <-
  hzar.copyModelLabels(rt$ror2$analysis$initDGs,
                       rt$ror2$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$ror2$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$ror2$runs$chains,
                                hzar.dataGroup.add),
                         rt$ror2$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$ror2$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$ror2$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$ror2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$ror2$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$ror2$analysis$model.name <-
        rownames(rt$ror2$analysis$AICcTable
        )[[ which.min(rt$ror2$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$ror2$analysis$model.selected <-
  rt$ror2$analysis$oDG$data.groups[[rt$ror2$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$ror2$analysis$model.selected,
                         names(rt$ror2$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$ror2$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$ror2$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$ror2$analysis$model.selected);

## End Analysis

#dev.off()

par(mfrow=c(4,1))
hzar.plot.fzCline(rt$back$analysis$model.selected);
hzar.plot.fzCline(rt$kitlg$analysis$model.selected);
hzar.plot.fzCline(rt$bnc2$analysis$model.selected);
hzar.plot.fzCline(rt$apc.camk4$analysis$model.selected);

### Analysis on APC/CAMK4 region (rustica-tytleri)--------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$apc.camk4 <- list();
## Space to hold the observed data
rt$apc.camk4$obs <- list();
## Space to hold the models to fit
rt$apc.camk4$models <- list();
## Space to hold the compiled fit requests
rt$apc.camk4$fitRs <- list();
## Space to hold the output data chains
rt$apc.camk4$runs <- list();
## Space to hold the analysed data
rt$apc.camk4$analysis <- list();

## Assign data for APC/CAMK4 region
rt$apc.camk4$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_apc_camk4_mean,
                             rt.dists$hi_apc_camk4_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$apc.camk4$obs);

## Make a helper function to set cline models
rt.loadapc.camk4model <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
  rt$apc.camk4$models[[id]] <<- hzar.makeCline1DFreq(rt$apc.camk4$obs, scaling, tails)

rt.loadapc.camk4model("fixed","none","modelI");
rt.loadapc.camk4model("free" ,"none","modelII");
rt.loadapc.camk4model("free" ,"both","modelIII");
rt.loadapc.camk4model("free" ,"right","modelIV");
rt.loadapc.camk4model("free" ,"left","modelV");
rt.loadapc.camk4model("free" ,"mirror","modelVI");

## Check the default settings
print(rt$apc.camk4$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$apc.camk4$models <- sapply(rt$apc.camk4$models,
                              hzar.model.addBoxReq,
                              -30 , 1973,
                              simplify=FALSE)

## Check the updated settings
print(rt$apc.camk4$models)

## Compile each of the models to prepare for fitting
rt$apc.camk4$fitRs$init <- sapply(rt$apc.camk4$models,
                                  hzar.first.fitRequest.old.ML,
                                  obsData=rt$apc.camk4$obs,
                                  verbose=FALSE,
                                  simplify=FALSE)

## Update the settings for the fitter if desired.
rt$apc.camk4$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$apc.camk4$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$apc.camk4$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$apc.camk4$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$apc.camk4$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$apc.camk4$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$apc.camk4$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$apc.camk4$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$apc.camk4$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$apc.camk4$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$apc.camk4$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$apc.camk4$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$apc.camk4$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$apc.camk4$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$apc.camk4$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$apc.camk4$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$apc.camk4$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$apc.camk4$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$apc.camk4$fitRs$init)

## Run the models for an initial chain
rt$apc.camk4$runs$init <- list()

rt$apc.camk4$runs$init$modelI <-
  hzar.doFit(rt$apc.camk4$fitRs$init$modelI)

rt$apc.camk4$runs$init$modelII <-
  hzar.doFit(rt$apc.camk4$fitRs$init$modelII)

rt$apc.camk4$runs$init$modelIII <-
  hzar.doFit(rt$apc.camk4$fitRs$init$modelIII)

rt$apc.camk4$runs$init$modelIV <-
  hzar.doFit(rt$apc.camk4$fitRs$init$modelIV)

rt$apc.camk4$runs$init$modelV <-
  hzar.doFit(rt$apc.camk4$fitRs$init$modelV)

rt$apc.camk4$runs$init$modelVI <-
  hzar.doFit(rt$apc.camk4$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$apc.camk4$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$apc.camk4$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$apc.camk4$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$apc.camk4$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$apc.camk4$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$apc.camk4$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$apc.camk4$fitRs$chains <-
  lapply(rt$apc.camk4$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$apc.camk4$fitRs$chains <-
  hzar.multiFitRequest(rt$apc.camk4$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$apc.camk4$runs$chains <-  hzar.doChain.multi(rt$apc.camk4$fitRs$chains,
                                                doPar=TRUE,
                                                inOrder=FALSE,
                                                count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$apc.camk4$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$apc.camk4$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$apc.camk4$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$apc.camk4$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$apc.camk4$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$apc.camk4$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$apc.camk4$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$apc.camk4$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$apc.camk4$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$apc.camk4$runs$init$modelI)
rt$apc.camk4$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$apc.camk4$runs$init$modelII)
rt$apc.camk4$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$apc.camk4$runs$init$modelIII)
rt$apc.camk4$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$apc.camk4$runs$init$modelIV)
rt$apc.camk4$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$apc.camk4$runs$init$modelV)
rt$apc.camk4$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$apc.camk4$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$apc.camk4$analysis$oDG <-
  hzar.make.obsDataGroup(rt$apc.camk4$analysis$initDGs)
rt$apc.camk4$analysis$oDG <-
  hzar.copyModelLabels(rt$apc.camk4$analysis$initDGs,
                       rt$apc.camk4$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$apc.camk4$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$apc.camk4$runs$chains,
                                hzar.dataGroup.add),
                         rt$apc.camk4$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$apc.camk4$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$apc.camk4$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$apc.camk4$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$apc.camk4$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$apc.camk4$analysis$model.name <-
        rownames(rt$apc.camk4$analysis$AICcTable
        )[[ which.min(rt$apc.camk4$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$apc.camk4$analysis$model.selected <-
  rt$apc.camk4$analysis$oDG$data.groups[[rt$apc.camk4$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$apc.camk4$analysis$model.selected,
                         names(rt$apc.camk4$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$apc.camk4$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$apc.camk4$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$apc.camk4$analysis$model.selected);

## End Analysis

#dev.off()

### Analysis on ICE1/lncRNA region (rustica-tytleri)--------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$ice1.lncrna <- list();
## Space to hold the observed data
rt$ice1.lncrna$obs <- list();
## Space to hold the models to fit
rt$ice1.lncrna$models <- list();
## Space to hold the compiled fit requests
rt$ice1.lncrna$fitRs <- list();
## Space to hold the output data chains
rt$ice1.lncrna$runs <- list();
## Space to hold the analysed data
rt$ice1.lncrna$analysis <- list();

## Assign data for APC/CAMK4 region
rt$ice1.lncrna$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_ice1_lncrna_mean,
                             rt.dists$hi_ice1_lncrna_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$ice1.lncrna$obs);

## Make a helper function to set cline models
rt.loadice1.lncrnamodel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
  rt$ice1.lncrna$models[[id]] <<- hzar.makeCline1DFreq(rt$ice1.lncrna$obs, scaling, tails)

rt.loadice1.lncrnamodel("fixed","none","modelI");
rt.loadice1.lncrnamodel("free" ,"none","modelII");
rt.loadice1.lncrnamodel("free" ,"both","modelIII");
rt.loadice1.lncrnamodel("free" ,"right","modelIV");
rt.loadice1.lncrnamodel("free" ,"left","modelV");
rt.loadice1.lncrnamodel("free" ,"mirror","modelVI");

## Check the default settings
print(rt$ice1.lncrna$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$ice1.lncrna$models <- sapply(rt$ice1.lncrna$models,
                              hzar.model.addBoxReq,
                              -30 , 1973,
                              simplify=FALSE)

## Check the updated settings
print(rt$ice1.lncrna$models)

## Compile each of the models to prepare for fitting
rt$ice1.lncrna$fitRs$init <- sapply(rt$ice1.lncrna$models,
                                  hzar.first.fitRequest.old.ML,
                                  obsData=rt$ice1.lncrna$obs,
                                  verbose=FALSE,
                                  simplify=FALSE)

## Update the settings for the fitter if desired.
rt$ice1.lncrna$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ice1.lncrna$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ice1.lncrna$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$ice1.lncrna$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ice1.lncrna$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ice1.lncrna$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$ice1.lncrna$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ice1.lncrna$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ice1.lncrna$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$ice1.lncrna$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ice1.lncrna$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ice1.lncrna$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$ice1.lncrna$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ice1.lncrna$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ice1.lncrna$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$ice1.lncrna$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$ice1.lncrna$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$ice1.lncrna$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$ice1.lncrna$fitRs$init)

## Run the models for an initial chain
rt$ice1.lncrna$runs$init <- list()

rt$ice1.lncrna$runs$init$modelI <-
  hzar.doFit(rt$ice1.lncrna$fitRs$init$modelI)

rt$ice1.lncrna$runs$init$modelII <-
  hzar.doFit(rt$ice1.lncrna$fitRs$init$modelII)

rt$ice1.lncrna$runs$init$modelIII <-
  hzar.doFit(rt$ice1.lncrna$fitRs$init$modelIII)

rt$ice1.lncrna$runs$init$modelIV <-
  hzar.doFit(rt$ice1.lncrna$fitRs$init$modelIV)

rt$ice1.lncrna$runs$init$modelV <-
  hzar.doFit(rt$ice1.lncrna$fitRs$init$modelV)

rt$ice1.lncrna$runs$init$modelVI <-
  hzar.doFit(rt$ice1.lncrna$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$ice1.lncrna$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$ice1.lncrna$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$ice1.lncrna$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$ice1.lncrna$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$ice1.lncrna$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$ice1.lncrna$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$ice1.lncrna$fitRs$chains <-
  lapply(rt$ice1.lncrna$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$ice1.lncrna$fitRs$chains <-
  hzar.multiFitRequest(rt$ice1.lncrna$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$ice1.lncrna$runs$chains <-  hzar.doChain.multi(rt$ice1.lncrna$fitRs$chains,
                                                doPar=TRUE,
                                                inOrder=FALSE,
                                                count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$ice1.lncrna$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$ice1.lncrna$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$ice1.lncrna$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$ice1.lncrna$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$ice1.lncrna$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$ice1.lncrna$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$ice1.lncrna$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$ice1.lncrna$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$ice1.lncrna$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$ice1.lncrna$runs$init$modelI)
rt$ice1.lncrna$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$ice1.lncrna$runs$init$modelII)
rt$ice1.lncrna$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$ice1.lncrna$runs$init$modelIII)
rt$ice1.lncrna$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$ice1.lncrna$runs$init$modelIV)
rt$ice1.lncrna$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$ice1.lncrna$runs$init$modelV)
rt$ice1.lncrna$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$ice1.lncrna$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$ice1.lncrna$analysis$oDG <-
  hzar.make.obsDataGroup(rt$ice1.lncrna$analysis$initDGs)
rt$ice1.lncrna$analysis$oDG <-
  hzar.copyModelLabels(rt$ice1.lncrna$analysis$initDGs,
                       rt$ice1.lncrna$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$ice1.lncrna$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$ice1.lncrna$runs$chains,
                                hzar.dataGroup.add),
                         rt$ice1.lncrna$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$ice1.lncrna$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$ice1.lncrna$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$ice1.lncrna$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$ice1.lncrna$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$ice1.lncrna$analysis$model.name <-
        rownames(rt$ice1.lncrna$analysis$AICcTable
        )[[ which.min(rt$ice1.lncrna$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$ice1.lncrna$analysis$model.selected <-
  rt$ice1.lncrna$analysis$oDG$data.groups[[rt$ice1.lncrna$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$ice1.lncrna$analysis$model.selected,
                         names(rt$ice1.lncrna$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$ice1.lncrna$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$ice1.lncrna$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$ice1.lncrna$analysis$model.selected);

## End Analysis

#dev.off()


### Analysis on PDE1C region (rustica-tytleri)--------------------------

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt$pde1c <- list();
## Space to hold the observed data
rt$pde1c$obs <- list();
## Space to hold the models to fit
rt$pde1c$models <- list();
## Space to hold the compiled fit requests
rt$pde1c$fitRs <- list();
## Space to hold the output data chains
rt$pde1c$runs <- list();
## Space to hold the analysed data
rt$pde1c$analysis <- list();

## Assign data for APC/CAMK4 region
rt$pde1c$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$hi_pde1c_mean,
                             rt.dists$hi_pde1c_count);

## Look at a graph of the observed data
hzar.plot.obsData(rt$pde1c$obs);

## Make a helper function to set cline models
rt.loadpde1cmodel <- function(scaling,tails,
                                  id=paste(scaling,tails,sep="."))
  rt$pde1c$models[[id]] <<- hzar.makeCline1DFreq(rt$pde1c$obs, scaling, tails)

rt.loadpde1cmodel("fixed","none","modelI");
rt.loadpde1cmodel("free" ,"none","modelII");
rt.loadpde1cmodel("free" ,"both","modelIII");
rt.loadpde1cmodel("free" ,"right","modelIV");
rt.loadpde1cmodel("free" ,"left","modelV");
rt.loadpde1cmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rt$pde1c$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt$pde1c$models <- sapply(rt$pde1c$models,
                              hzar.model.addBoxReq,
                              -30 , 1973,
                              simplify=FALSE)

## Check the updated settings
print(rt$pde1c$models)

## Compile each of the models to prepare for fitting
rt$pde1c$fitRs$init <- sapply(rt$pde1c$models,
                                  hzar.first.fitRequest.old.ML,
                                  obsData=rt$pde1c$obs,
                                  verbose=FALSE,
                                  simplify=FALSE)

## Update the settings for the fitter if desired.
rt$pde1c$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$pde1c$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$pde1c$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt$pde1c$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$pde1c$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$pde1c$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt$pde1c$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$pde1c$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$pde1c$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$pde1c$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$pde1c$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$pde1c$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$pde1c$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$pde1c$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$pde1c$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt$pde1c$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt$pde1c$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt$pde1c$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt$pde1c$fitRs$init)

## Run the models for an initial chain
rt$pde1c$runs$init <- list()

rt$pde1c$runs$init$modelI <-
  hzar.doFit(rt$pde1c$fitRs$init$modelI)

rt$pde1c$runs$init$modelII <-
  hzar.doFit(rt$pde1c$fitRs$init$modelII)

rt$pde1c$runs$init$modelIII <-
  hzar.doFit(rt$pde1c$fitRs$init$modelIII)

rt$pde1c$runs$init$modelIV <-
  hzar.doFit(rt$pde1c$fitRs$init$modelIV)

rt$pde1c$runs$init$modelV <-
  hzar.doFit(rt$pde1c$fitRs$init$modelV)

rt$pde1c$runs$init$modelVI <-
  hzar.doFit(rt$pde1c$fitRs$init$modelVI)

## Un-comment to plot the traces
#plot(hzar.mcmc.bindLL(rt$pde1c$runs$init$modelI))
#plot(hzar.mcmc.bindLL(rt$pde1c$runs$init$modelII))
#plot(hzar.mcmc.bindLL(rt$pde1c$runs$init$modelIII))
#plot(hzar.mcmc.bindLL(rt$pde1c$runs$init$modelIV))
#plot(hzar.mcmc.bindLL(rt$pde1c$runs$init$modelV))
#plot(hzar.mcmc.bindLL(rt$pde1c$runs$init$modelVI))

## Compile a new set of fit requests using the initial chains 
rt$pde1c$fitRs$chains <-
  lapply(rt$pde1c$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt$pde1c$fitRs$chains <-
  hzar.multiFitRequest(rt$pde1c$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt$pde1c$runs$chains <-  hzar.doChain.multi(rt$pde1c$fitRs$chains,
                                                doPar=TRUE,
                                                inOrder=FALSE,
                                                count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt$pde1c$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt$pde1c$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt$pde1c$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt$pde1c$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt$pde1c$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt$pde1c$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt$pde1c$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt$pde1c$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt$pde1c$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt$pde1c$runs$init$modelI)
rt$pde1c$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt$pde1c$runs$init$modelII)
rt$pde1c$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt$pde1c$runs$init$modelIII)
rt$pde1c$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt$pde1c$runs$init$modelIV)
rt$pde1c$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt$pde1c$runs$init$modelV)
rt$pde1c$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt$pde1c$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt$pde1c$analysis$oDG <-
  hzar.make.obsDataGroup(rt$pde1c$analysis$initDGs)
rt$pde1c$analysis$oDG <-
  hzar.copyModelLabels(rt$pde1c$analysis$initDGs,
                       rt$pde1c$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt$pde1c$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt$pde1c$runs$chains,
                                hzar.dataGroup.add),
                         rt$pde1c$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt$pde1c$analysis$oDG$data.groups))

## Compare the cline models to the null model graphically
hzar.plot.cline(rt$pde1c$analysis$oDG);

## Do model selection based on the AICc scores
print(rt$pde1c$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$pde1c$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt$pde1c$analysis$model.name <-
        rownames(rt$pde1c$analysis$AICcTable
        )[[ which.min(rt$pde1c$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt$pde1c$analysis$model.selected <-
  rt$pde1c$analysis$oDG$data.groups[[rt$pde1c$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$pde1c$analysis$model.selected,
                         names(rt$pde1c$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$pde1c$analysis$model.selected))

## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$pde1c$analysis$model.selected);

hzar.plot.cline(rt$vent$analysis$model.selected);

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(rt$pde1c$analysis$model.selected);

## End Analysis

#dev.off()

### Save analysis results to Rdata object-----------------------------------

save(rt, file = './candidate_results_hzar/hzar_data_rustica-tytleri.RData')
