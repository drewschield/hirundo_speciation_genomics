#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Note: this script expects a command line argument in the form of 'chr1A', for example.

## Load libraries

library(tidyverse)
library(data.table)
library(R.utils)
library(vcfR)
library(rehh)

## Convert VCF to haplo format

ru.chrom <- data2haplohh(hap_file = paste0('./shapeit/results/hirundo_rustica.parental.rustica.',args[1],'.snps.miss05.maf01.phased.vcf.gz'), polarize_vcf = FALSE)
ty.chrom <- data2haplohh(hap_file = paste0('./shapeit/results/hirundo_rustica.parental.tytleri.',args[1],'.snps.miss05.maf01.phased.vcf.gz'), polarize_vcf = FALSE)
gu.chrom <- data2haplohh(hap_file = paste0('./shapeit/results/hirundo_rustica.parental.gutturalis.',args[1],'.snps.miss05.maf01.phased.vcf.gz'), polarize_vcf = FALSE)

## Filter SNPs with MAF < 0.05

ru.chrom_f <- subset(ru.chrom, min_maf = 0.05)
ty.chrom_f <- subset(ty.chrom, min_maf = 0.05)
gu.chrom_f <- subset(gu.chrom, min_maf = 0.05)

## Run scans

ru.chrom_scan <- scan_hh(ru.chrom_f, polarized = FALSE)
ty.chrom_scan <- scan_hh(ty.chrom_f, polarized = FALSE)
gu.chrom_scan <- scan_hh(gu.chrom_f, polarized = FALSE)

## Calculate iHS

ru.chrom_ihs <- ihh2ihs(ru.chrom_scan, freqbin = 1)
ty.chrom_ihs <- ihh2ihs(ty.chrom_scan, freqbin = 1)
gu.chrom_ihs <- ihh2ihs(gu.chrom_scan, freqbin = 1)

## Calculate xp-EHH

ruty.chrom_xpehh <- ies2xpehh(ru.chrom_scan, ty.chrom_scan, popname1 = "rustica", popname2 = "tytleri", include_freq = T)
rugu.chrom_xpehh <- ies2xpehh(ru.chrom_scan, gu.chrom_scan, popname1 = "rustica", popname2 = "gutturalis", include_freq = T)
guty.chrom_xpehh <- ies2xpehh(ty.chrom_scan, gu.chrom_scan, popname1 = "tytleri", popname2 = "gutturalis", include_freq = T)

## Save scans to Rdata objects

save(ru.chrom_scan, file = paste0('./Rdata/rehh_scan_',args[1],'_rustica.RData'))
save(ty.chrom_scan, file = paste0('./Rdata/rehh_scan_',args[1],'_tytleri.RData'))
save(gu.chrom_scan, file = paste0('./Rdata/rehh_scan_',args[1],'_gutturalis.RData'))

## Write iHS results to files

write.table(ru.chrom_ihs$ihs,file=paste0('./results/iHS_',args[1],'_rustica.txt'),row.names=F,quote=F,sep="\t")
write.table(ty.chrom_ihs$ihs,file=paste0('./results/iHS_',args[1],'_tytleri.txt'),row.names=F,quote=F,sep="\t")
write.table(gu.chrom_ihs$ihs,file=paste0('./results/iHS_',args[1],'_gutturalis.txt'),row.names=F,quote=F,sep="\t")

## Write xp-EHH results to files

write.table(ruty.chrom_xpehh,file=paste0('./results/XP-EHH_',args[1],'_rustica-tytleri.txt'),row.names=F,quote=F,sep="\t")
write.table(rugu.chrom_xpehh,file=paste0('./results/XP-EHH_',args[1],'_rustica-gutturalis.txt'),row.names=F,quote=F,sep="\t")
write.table(guty.chrom_xpehh,file=paste0('./results/XP-EHH_',args[1],'_tytleri-gutturalis.txt'),row.names=F,quote=F,sep="\t")

## End Analysis

quit(save="no")

