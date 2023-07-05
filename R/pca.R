############################################################################
# Principal components analysis of Barn Swallow SNPs
############################################################################

### Goal: compare PCA clustering of individuals based on autosomal & Z-linked
### SNPs.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory---------------------------------------------------

setwd('~/hirundo_speciation_genomics/analysis/pca/')

### Install SNPRelate and other dependencies--------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
library(SNPRelate)
library(scales)
library(RColorBrewer)

### ------------------------------------------------------------------------
### Analysis of autosomal + Z-linked SNPs
### ------------------------------------------------------------------------

genofile_hirundoA <- snpgdsVCF2GDS('./hirundo_rustica+smithii.allsites.final.snps.miss02.maf1.ingroup.indv.auto+chrZ.vcf.gz','./hirundo_rustica.auto+chrZ.gds')
genofile_hirundoA <- openfn.gds('hirundo_rustica.auto+chrZ.gds')
hirundoA_pca <- snpgdsPCA(genofile_hirundoA,autosome.only = FALSE)
write.table(hirundoA_pca$eigenvect,file = "hirundo_rustica.auto+chrZ.pca_eigenvect.txt",quote=FALSE,row.names = FALSE)

### [NOTE]:
### Append relevant metadata to the tab-delimited table of PC loadings and read the file back in...
### Read table of eigenvectors (not directly relevant unless you've appended the sample metadata to the table of eigenvectors)

hirundo_chromA <- read.table("hirundo_rustica.auto+chrZ.pca_eigenvect+metadata.txt",header=T)

### Calculate variance------------------------------------------------------

pc.percent_hirundoA <- hirundoA_pca$varprop*100
head(round(pc.percent_hirundoA, 2))

plot(pc.percent_hirundoA)

### Set PCA objects---------------------------------------------------------

PC1_hrA <- hirundo_chromA$V1
PC2_hrA <- hirundo_chromA$V2
ID_hrA <- hirundo_chromA$sample
pop_hrA <- hirundo_chromA$pop

### Plot PCs----------------------------------------------------------------

palette("default")
cc <- palette()
palette(c(cc,"purple","brown"))
palette()
par(mfrow=c(1,1))
plot(PC1_hrA,PC2_hrA,pch=20,col=hirundo_chromA$pop,xlab='PC1 (2.85%)',ylab='PC2 (1.33%)',main='Autosomes + Z chromosome')
legend('topleft', legend = levels(hirundo_chromZ$pop), col = 1:10, cex = 0.8, pch = 20)

#plot with reverse x-axis to match map
plot(PC1_hrA,PC2_hrA,pch=20,col=hirundo_chromA$pop,xlab='PC1 (2.85%)',ylab='PC2 (1.33%)',main='Autosomes + Z chromosome',xlim=c(0.1,-0.1))
