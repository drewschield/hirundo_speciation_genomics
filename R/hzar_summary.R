############################################################################
# Barn swallow geographic clines in hybrid zone transects
############################################################################

# This script contains commands to summarize the results of geographic
# cline analyses in HZAR. Here, we will look at model AICc support and 
# parameters for comparison and for reporting the results.

# This workflow involves:
# 1. Reading in and summarizing genome background cline parameters; outputting results
# 2. Summarizing of AICc tables and cline center/width summary statistics for candidate loci

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/hzar/')

#install.packages("hzar") --> had to install from source because HZAR has been removed from CRAN(?)
#install.packages('./hzar_0.2-5.tar.gz', repos = NULL, type ="source")

library(Rmisc)
library(tidyverse)
library(data.table)
library(hzar)
library(doMC)

### Load data for candidate loci-------------------------------------------

load("./candidate_results_hzar/hzar_data_rustica-tytleri.RData")
load("./candidate_results_hzar/hzar_data_rustica-gutturalis.RData")
load("./candidate_results_hzar/hzar_data_tytleri-gutturalis.RData")

## Set as other variables (background loci have same list names)

rt.cand <- rt
rg.cand <- rg
tg.cand <- tg

### Read in background locus lists-----------------------------------------

rt.list <- read.table('list.background_locus.rustica-tytleri.txt',header=F)
rg.list <- read.table('list.background_locus.rustica-gutturalis.txt',header=F)
tg.list <- read.table('list.background_locus.tytleri-gutturalis.txt',header=F)

### Extract genome background cline parameters------------------------------

# rustica-tytleri
## Make empty vectors
rt.back.center <- c()
rt.back.width <- c()
rt.back.pmin <- c()
rt.back.pmax <- c()
rt.back.deltah <- c()
## Add data per locus to vectors
for (i in rt.list$V1){
  print(i)
  load(file=paste0('./results_hzar/hzar_rustica-tytleri.background_',i,'.RData'))
  center <- rt[[1]]$analysis$model.selected$ML.cline$param.all$center
  width <- rt[[1]]$analysis$model.selected$ML.cline$param.all$width
  pmin <- rt[[1]]$analysis$model.selected$ML.cline$param.all$pMin
  pmax <- rt[[1]]$analysis$model.selected$ML.cline$param.all$pMax
  deltah <- rt[[1]]$analysis$model.selected$ML.cline$param.all$pMax - rt[[1]]$analysis$model.selected$ML.cline$param.all$pMin
  rt.back.center=append(rt.back.center,center)
  rt.back.width=append(rt.back.width,width)
  rt.back.pmin=append(rt.back.pmin,pmin)
  rt.back.pmax=append(rt.back.pmax,pmax)
  rt.back.deltah=append(rt.back.deltah,deltah)
}
## Make empty dataframe and add to it
rt.back <- c()
rt.back$center <- rt.back.center
rt.back$width <- rt.back.width
rt.back$pmin <- rt.back.pmin
rt.back$pmax <- rt.back.pmax
rt.back$deltah <- rt.back.deltah
write.table(rt.back,file='./background_locus.param.rustica-tytleri.txt',quote=F,row.names=F,sep="\t")

# rustica-gutturalis
## Make empty vectors
rg.back.center <- c()
rg.back.width <- c()
rg.back.pmin <- c()
rg.back.pmax <- c()
rg.back.deltah <- c()
## Add data per locus to vectors
for (i in rg.list$V1){
  print(i)
  load(file=paste0('./results_hzar/hzar_rustica-gutturalis.background_',i,'.RData'))
  center <- rg[[1]]$analysis$model.selected$ML.cline$param.all$center
  width <- rg[[1]]$analysis$model.selected$ML.cline$param.all$width
  pmin <- rg[[1]]$analysis$model.selected$ML.cline$param.all$pMin
  pmax <- rg[[1]]$analysis$model.selected$ML.cline$param.all$pMax
  deltah <- rg[[1]]$analysis$model.selected$ML.cline$param.all$pMax - rg[[1]]$analysis$model.selected$ML.cline$param.all$pMin
  rg.back.center=append(rg.back.center,center)
  rg.back.width=append(rg.back.width,width)
  rg.back.pmin=append(rg.back.pmin,pmin)
  rg.back.pmax=append(rg.back.pmax,pmax)
  rg.back.deltah=append(rg.back.deltah,deltah)
}
## Make empty dataframe and add to it
rg.back <- c()
rg.back$center <- rg.back.center
rg.back$width <- rg.back.width
rg.back$pmin <- rg.back.pmin
rg.back$pmax <- rg.back.pmax
rg.back$deltah <- rg.back.deltah
write.table(rg.back,file='./background_locus.param.rustica-gutturalis.txt',quote=F,row.names=F,sep="\t")

# tytleri-gutturalis
## Make empty vectors
tg.back.center <- c()
tg.back.width <- c()
tg.back.pmin <- c()
tg.back.pmax <- c()
tg.back.deltah <- c()
## Add data per locus to vectors
for (i in tg.list$V1){
  print(i)
  load(file=paste0('./results_hzar/hzar_tytleri-gutturalis.background_',i,'.RData'))
  center <- tg[[1]]$analysis$model.selected$ML.cline$param.all$center
  width <- tg[[1]]$analysis$model.selected$ML.cline$param.all$width
  pmin <- tg[[1]]$analysis$model.selected$ML.cline$param.all$pMin
  pmax <- tg[[1]]$analysis$model.selected$ML.cline$param.all$pMax
  deltah <- tg[[1]]$analysis$model.selected$ML.cline$param.all$pMax - tg[[1]]$analysis$model.selected$ML.cline$param.all$pMin
  tg.back.center=append(tg.back.center,center)
  tg.back.width=append(tg.back.width,width)
  tg.back.pmin=append(tg.back.pmin,pmin)
  tg.back.pmax=append(tg.back.pmax,pmax)
  tg.back.deltah=append(tg.back.deltah,deltah)
}
## Make empty dataframe and add to it
tg.back <- c()
tg.back$center <- tg.back.center
tg.back$width <- tg.back.width
tg.back$pmin <- tg.back.pmin
tg.back$pmax <- tg.back.pmax
tg.back$deltah <- tg.back.deltah
write.table(tg.back,file='./background_locus.param.tytleri-gutturalis.txt',quote=F,row.names=F,sep="\t")

### Read in genome background parameters------------------------------------

## If above has already ran!
rt.back <- read.table('./background_locus.param.rustica-tytleri.txt',header=T)
rg.back <- read.table('./background_locus.param.rustica-gutturalis.txt',header=T)
tg.back <- read.table('./background_locus.param.tytleri-gutturalis.txt',header=T)

### INSERT CANDIDATE LOCUS MODEL SELECTION HERE----

### SUMMARY OF BACKGROUND-----

mean(rt.back$center)
sd(rt.back$center)
mean(rt.back$width)
sd(rt.back$width)
mean(rt.back$pmin)
sd(rt.back$pmin)
mean(rt.back$pmax)
sd(rt.back$pmax)

mean(rg.back$center)
sd(rg.back$center)
mean(rg.back$width)
sd(rg.back$width)
mean(rg.back$pmin)
sd(rg.back$pmin)
mean(rg.back$pmax)
sd(rg.back$pmax)

mean(tg.back$center)
sd(tg.back$center)
mean(tg.back$width)
sd(tg.back$width)
mean(tg.back$pmin)
sd(tg.back$pmin)
mean(tg.back$pmax)
sd(tg.back$pmax)

### Cline model parameters for candidate loci--------------------------------

# rustica-tytleri
## First, set up empty vectors for variables, then loop through selected model parameters
loci <- c('kitlg','plxnc1','spef2.prlr','slc45a2','bnc2','gnaq','ror2','apc.camk4','ice1.lncrna','pde1c')
rt.loci.locus <- c()
rt.loci.center <- c()
rt.loci.width <- c()
rt.loci.pmin <- c()
rt.loci.pmax <- c()
for (locus in loci){
  rt.loci.locus=append(rt.loci.locus,locus)
  rt.loci.center=append(rt.loci.center, rt[[locus]]$analysis$model.selected$ML.cline$param.all$center)
  rt.loci.width=append(rt.loci.width, rt[[locus]]$analysis$model.selected$ML.cline$param.all$width)
  rt.loci.pmin=append(rt.loci.pmin, rt[[locus]]$analysis$model.selected$ML.cline$param.all$pMin)
  rt.loci.pmax=append(rt.loci.pmax, rt[[locus]]$analysis$model.selected$ML.cline$param.all$pMax)
}

## Some of the loci had similar support for model II, which was a much better visual fit to the data. Replace with these values
rt.loci.center[3]=replace(rt.loci.center[3],values=rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rt.loci.width[3]=replace(rt.loci.width[3],values=rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rt.loci.pmin[3]=replace(rt.loci.pmin[3],values=rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rt.loci.pmax[3]=replace(rt.loci.pmax[3],values=rt$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

rt.loci.center[4]=replace(rt.loci.center[4],values=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rt.loci.width[4]=replace(rt.loci.width[4],values=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rt.loci.pmin[4]=replace(rt.loci.pmin[4],values=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rt.loci.pmax[4]=replace(rt.loci.pmax[4],values=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

rt.loci.center[7]=replace(rt.loci.center[7],values=rt$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rt.loci.width[7]=replace(rt.loci.width[7],values=rt$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rt.loci.pmin[7]=replace(rt.loci.pmin[7],values=rt$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rt.loci.pmax[7]=replace(rt.loci.pmax[7],values=rt$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

rt.loci.center[8]=replace(rt.loci.center[8],values=rt$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rt.loci.width[8]=replace(rt.loci.width[8],values=rt$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rt.loci.pmin[8]=replace(rt.loci.pmin[8],values=rt$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rt.loci.pmax[8]=replace(rt.loci.pmax[8],values=rt$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

rt.loci.center[11]=replace(rt.loci.center[11],values=rt$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rt.loci.width[11]=replace(rt.loci.width[11],values=rt$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rt.loci.pmin[11]=replace(rt.loci.pmin[11],values=rt$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rt.loci.pmax[11]=replace(rt.loci.pmax[11],values=rt$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

## Write parameters to file
rt.loci <- c()
rt.loci$locus <- rt.loci.locus
rt.loci$center <- rt.loci.center
rt.loci$width <- rt.loci.width
rt.loci$pmin <- rt.loci.pmin
rt.loci$pmax <- rt.loci.pmax
write.table(rt.loci,file='./candidate_locus.param.rustica-tytleri.txt',quote=F,row.names=F,sep="\t")

## Now gimme the other model parameters
loci <- c('spef2.prlr','slc45a2','ror2','apc.camk4','pde1c')
rt.loci.locus <- c()
rt.loci.center <- c()
rt.loci.width <- c()
rt.loci.pmin <- c()
rt.loci.pmax <- c()
for (locus in loci){
  rt.loci.locus=append(rt.loci.locus,locus)
  rt.loci.center=append(rt.loci.center, rt[[locus]]$analysis$model.selected$ML.cline$param.all$center)
  rt.loci.width=append(rt.loci.width, rt[[locus]]$analysis$model.selected$ML.cline$param.all$width)
  rt.loci.pmin=append(rt.loci.pmin, rt[[locus]]$analysis$model.selected$ML.cline$param.all$pMin)
  rt.loci.pmax=append(rt.loci.pmax, rt[[locus]]$analysis$model.selected$ML.cline$param.all$pMax)
}

rt.loci <- c()
rt.loci$locus <- rt.loci.locus
rt.loci$center <- rt.loci.center
rt.loci$width <- rt.loci.width
rt.loci$pmin <- rt.loci.pmin
rt.loci$pmax <- rt.loci.pmax
write.table(rt.loci,file='./candidate_locus.param.rustica-tytleri_other_model.txt',quote=F,row.names=F,sep="\t")

# rustica-gutturalis
## First, set up empty vectors for variables, then loop through selected model parameters
loci <- c('kitlg','plxnc1','spef2.prlr','slc45a2','bnc2','gnaq','ror2','apc.camk4','ice1.lncrna','pde1c')
rg.loci.locus <- c()
rg.loci.center <- c()
rg.loci.width <- c()
rg.loci.pmin <- c()
rg.loci.pmax <- c()
for (locus in loci){
  rg.loci.locus=append(rg.loci.locus,locus)
  rg.loci.center=append(rg.loci.center, rg[[locus]]$analysis$model.selected$ML.cline$param.all$center)
  rg.loci.width=append(rg.loci.width, rg[[locus]]$analysis$model.selected$ML.cline$param.all$width)
  rg.loci.pmin=append(rg.loci.pmin, rg[[locus]]$analysis$model.selected$ML.cline$param.all$pMin)
  rg.loci.pmax=append(rg.loci.pmax, rg[[locus]]$analysis$model.selected$ML.cline$param.all$pMax)
}

## Some of the loci had similar supporg for model II, which was a much better visual fit to the data. Replace with these values
rg.loci.center[4]=replace(rg.loci.center[4],values=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rg.loci.width[4]=replace(rg.loci.width[4],values=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rg.loci.pmin[4]=replace(rg.loci.pmin[4],values=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rg.loci.pmax[4]=replace(rg.loci.pmax[4],values=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

rg.loci.center[9]=replace(rg.loci.center[9],values=rg$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$center)
rg.loci.width[9]=replace(rg.loci.width[9],values=rg$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)
rg.loci.pmin[9]=replace(rg.loci.pmin[9],values=rg$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin)
rg.loci.pmax[9]=replace(rg.loci.pmax[9],values=rg$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax)

## Write parameters to file
rg.loci <- c()
rg.loci$locus <- rg.loci.locus
rg.loci$center <- rg.loci.center
rg.loci$width <- rg.loci.width
rg.loci$pmin <- rg.loci.pmin
rg.loci$pmax <- rg.loci.pmax
write.table(rg.loci,file='./candidate_locus.param.rustica-gutturalis.txt',quote=F,row.names=F,sep="\t")

## Now gimme the other model parameters
loci <- c('slc45a2','ice1.lncrna')
rg.loci.locus <- c()
rg.loci.center <- c()
rg.loci.width <- c()
rg.loci.pmin <- c()
rg.loci.pmax <- c()
for (locus in loci){
  rg.loci.locus=append(rg.loci.locus,locus)
  rg.loci.center=append(rg.loci.center, rg[[locus]]$analysis$model.selected$ML.cline$param.all$center)
  rg.loci.width=append(rg.loci.width, rg[[locus]]$analysis$model.selected$ML.cline$param.all$width)
  rg.loci.pmin=append(rg.loci.pmin, rg[[locus]]$analysis$model.selected$ML.cline$param.all$pMin)
  rg.loci.pmax=append(rg.loci.pmax, rg[[locus]]$analysis$model.selected$ML.cline$param.all$pMax)
}

rg.loci <- c()
rg.loci$locus <- rg.loci.locus
rg.loci$center <- rg.loci.center
rg.loci$width <- rg.loci.width
rg.loci$pmin <- rg.loci.pmin
rg.loci$pmax <- rg.loci.pmax
write.table(rg.loci,file='./candidate_locus.param.rustica-gutturalis_other_model.txt',quote=F,row.names=F,sep="\t")

# tytleri-gutturalis
## First, set up empty vectors for variables, then loop through selected model parameters
loci <- c('kitlg','plxnc1','spef2.prlr','slc45a2','bnc2','gnaq','ror2','apc.camk4','ice1.lncrna','pde1c')
tg.loci.locus <- c()
tg.loci.center <- c()
tg.loci.width <- c()
tg.loci.pmin <- c()
tg.loci.pmax <- c()
for (locus in loci){
  tg.loci.locus=append(tg.loci.locus,locus)
  tg.loci.center=append(tg.loci.center, tg[[locus]]$analysis$model.selected$ML.cline$param.all$center)
  tg.loci.width=append(tg.loci.width, tg[[locus]]$analysis$model.selected$ML.cline$param.all$width)
  tg.loci.pmin=append(tg.loci.pmin, tg[[locus]]$analysis$model.selected$ML.cline$param.all$pMin)
  tg.loci.pmax=append(tg.loci.pmax, tg[[locus]]$analysis$model.selected$ML.cline$param.all$pMax)
}

## Write parameters to file
tg.loci <- c()
tg.loci$locus <- tg.loci.locus
tg.loci$center <- tg.loci.center
tg.loci$width <- tg.loci.width
tg.loci$pmin <- tg.loci.pmin
tg.loci$pmax <- tg.loci.pmax
write.table(tg.loci,file='./candidate_locus.param.tytleri-gutturalis.txt',quote=F,row.names=F,sep="\t")


### Statistical comparison between background and candidates-----------------

# rustica-tytleri
## Cline center (c)
test <- c()
test=append(test,rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$param.all$center)
test=append(test,rt.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$center)
test=append(test,rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rt.cand$ror2$analysis$oDG$data.groups$modelII$param.all$center)
test=append(test,rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$param.all$center)
test=append(test,rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rt.cand$pde1c$analysis$oDG$data.groups$modelII$param.all$center)

wilcox.test(test,rt.back$center)

## Cline width (w)
test <- c()
test=append(test,rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$param.all$width)
test=append(test,rt.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$width)
test=append(test,rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rt.cand$ror2$analysis$oDG$data.groups$modelII$param.all$width)
test=append(test,rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$param.all$width)
test=append(test,rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rt.cand$pde1c$analysis$oDG$data.groups$modelII$param.all$width)

wilcox.test(test,rt.back$width)

## Cline delta h
test <- c()
test=append(test,rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMax - rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMax - rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$param.all$pMax - rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$param.all$pMin)
test=append(test,rt.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$pMax - rt.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$pMin)
test=append(test,rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMax - rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMax - rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rt.cand$ror2$analysis$oDG$data.groups$modelII$param.all$pMax - rt.cand$ror2$analysis$oDG$data.groups$modelII$param.all$pMin)
test=append(test,rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$param.all$pMax - rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$param.all$pMin)
test=append(test,rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMax - rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rt.cand$pde1c$analysis$oDG$data.groups$modelII$param.all$pMax - rt.cand$pde1c$analysis$oDG$data.groups$modelII$param.all$pMin)

wilcox.test(test,rt.back$pmax-rt.back$pmin)

# rustica-gutturalis
## Cline center (c)
test <- c()
test=append(test,rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$center)
test=append(test,rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$ror2$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center)
test=append(test,rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$param.all$center)
test=append(test,rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center)

wilcox.test(test,rg.back$center)

## Cline width (w)
test <- c()
test=append(test,rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$width)
test=append(test,rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$ror2$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$width)
test=append(test,rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$param.all$width)
test=append(test,rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$width)

wilcox.test(test,rg.back$width)

## Cline delta h
test <- c()
test=append(test,rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$pMax - rg.cand$slc45a2$analysis$oDG$data.groups$modelII$param.all$pMin)
test=append(test,rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$param.all$pMax - rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$param.all$pMin)
test=append(test,rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMax - rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMin)

wilcox.test(test,rg.back$pmax-rg.back$pmin)


# tytleri-gutturalis
## Cline center (c)
test <- c()
test=append(test,tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$ror2$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center)
test=append(test,tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center)

wilcox.test(test,tg.back$center)

## Cline width (w)
test <- c()
test=append(test,tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$ror2$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$width)
test=append(test,tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$width)

wilcox.test(test,tg.back$width)

## Cline delta h
test <- c()
test=append(test,tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMin)
test=append(test,tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMax - tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMin)

wilcox.test(test,tg.back$pmax-tg.back$pmin)


### Summaries for rustica-tytleri clines-----------------------------------

## KITLG

# Do model selection based on the AICc scores
print(rt$kitlg$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$kitlg$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$kitlg$analysis$model.name <-
        rownames(rt$kitlg$analysis$AICcTable
        )[[ which.min(rt$kitlg$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$kitlg$analysis$model.selected <-
  rt$kitlg$analysis$oDG$data.groups[[rt$kitlg$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$kitlg$analysis$model.selected,
                         names(rt$kitlg$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$kitlg$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$kitlg$analysis$model.selected, add=T,xlim=c(0,3500))

## SPEF2-PRLR

# Do model selection based on the AICc scores
print(rt$spef2.prlr$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$spef2.prlr$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$spef2.prlr$analysis$model.name <-
        rownames(rt$spef2.prlr$analysis$AICcTable
        )[[ which.min(rt$spef2.prlr$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$spef2.prlr$analysis$model.selected <-
  rt$spef2.prlr$analysis$oDG$data.groups[[rt$spef2.prlr$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$spef2.prlr$analysis$model.selected,
                         names(rt$spef2.prlr$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
#print(hzar.get.ML.cline(rt$spef2.prlr$analysis$model.selected))
# Plot the maximum likelihood cline for the selected model
#hzar.plot.cline(rt$spef2.prlr$analysis$model.selected, add=T, lty=2,xlim=c(0,3500));

# Similar AICc support for Model I and II; look at model II
print(hzar.get.ML.cline(rt$spef2.prlr$analysis$oDG$data.groups$modelII))
hzar.plot.cline(rt$spef2.prlr$analysis$oDG$data.groups$modelII, add=T,xlim=c(0,3500))

## SLC45A2

# Do model selection based on the AICc scores
print(rt$slc45a2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$slc45a2$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$slc45a2$analysis$model.name <-
        rownames(rt$slc45a2$analysis$AICcTable
        )[[ which.min(rt$slc45a2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$slc45a2$analysis$model.selected <-
  rt$slc45a2$analysis$oDG$data.groups[[rt$slc45a2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$slc45a2$analysis$model.selected,
                         names(rt$slc45a2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
#print(hzar.get.ML.cline(rt$slc45a2$analysis$model.selected))
# Plot the maximum likelihood cline for the selected model
#hzar.plot.cline(rt$slc45a2$analysis$model.selected, add=T, lty=2,xlim=c(0,3500));

# Similar AICc support for Model I and II; look at model II
print(hzar.get.ML.cline(rt$slc45a2$analysis$oDG$data.groups$modelII))
hzar.plot.cline(rt$slc45a2$analysis$oDG$data.groups$modelII, add=T,xlim=c(0,3500))

## BNC2

# Do model selection based on the AICc scores
print(rt$bnc2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$bnc2$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$bnc2$analysis$model.name <-
        rownames(rt$bnc2$analysis$AICcTable
        )[[ which.min(rt$bnc2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$bnc2$analysis$model.selected <-
  rt$bnc2$analysis$oDG$data.groups[[rt$bnc2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$bnc2$analysis$model.selected,
                         names(rt$bnc2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$bnc2$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$bnc2$analysis$model.selected, add=T,xlim=c(0,3500))

## ROR2

# Do model selection based on the AICc scores
print(rt$ror2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$ror2$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$ror2$analysis$model.name <-
        rownames(rt$ror2$analysis$AICcTable
        )[[ which.min(rt$ror2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$ror2$analysis$model.selected <-
  rt$ror2$analysis$oDG$data.groups[[rt$ror2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$ror2$analysis$model.selected,
                         names(rt$ror2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$ror2$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$ror2$analysis$model.selected, add=T,xlim=c(0,3500))

## APC-CAMK4

# Do model selection based on the AICc scores
print(rt$apc.camk4$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$apc.camk4$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$apc.camk4$analysis$model.name <-
        rownames(rt$apc.camk4$analysis$AICcTable
        )[[ which.min(rt$apc.camk4$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$apc.camk4$analysis$model.selected <-
  rt$apc.camk4$analysis$oDG$data.groups[[rt$apc.camk4$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$apc.camk4$analysis$model.selected,
                         names(rt$apc.camk4$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$apc.camk4$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$apc.camk4$analysis$model.selected, add=T,xlim=c(0,3500))

## PDE1C

# Do model selection based on the AICc scores
print(rt$pde1c$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt$pde1c$analysis$oDG));
# Print out the model with the minimum AICc score
print(rt$pde1c$analysis$model.name <-
        rownames(rt$pde1c$analysis$AICcTable
        )[[ which.min(rt$pde1c$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rt$pde1c$analysis$model.selected <-
  rt$pde1c$analysis$oDG$data.groups[[rt$pde1c$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt$pde1c$analysis$model.selected,
                         names(rt$pde1c$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt$pde1c$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rt$pde1c$analysis$model.selected, add=T,xlim=c(0,3500))


### Summaries for rustica-gutturalis clines--------------------------------

## KITLG

# Do model selection based on the AICc scores
print(rg$kitlg$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$kitlg$analysis$oDG));
# Print out the model with the minimum AICc score
print(rg$kitlg$analysis$model.name <-
        rownames(rg$kitlg$analysis$AICcTable
        )[[ which.min(rg$kitlg$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rg$kitlg$analysis$model.selected <-
  rg$kitlg$analysis$oDG$data.groups[[rg$kitlg$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$kitlg$analysis$model.selected,
                         names(rg$kitlg$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$kitlg$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$kitlg$analysis$model.selected, add=T,xlim=c(0,3500))

## SLC45A2

# Do model selection based on the AICc scores
print(rg$slc45a2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$slc45a2$analysis$oDG));
# Print out the model with the minimum AICc score
print(rg$slc45a2$analysis$model.name <-
        rownames(rg$slc45a2$analysis$AICcTable
        )[[ which.min(rg$slc45a2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rg$slc45a2$analysis$model.selected <-
  rg$slc45a2$analysis$oDG$data.groups[[rg$slc45a2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$slc45a2$analysis$model.selected,
                         names(rg$slc45a2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$slc45a2$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$slc45a2$analysis$model.selected, add=T, lty=2,xlim=c(0,3500));

# Similar AICc support for Model I and II; look at model II
print(hzar.get.ML.cline(rg$slc45a2$analysis$oDG$data.groups$modelII))
hzar.plot.cline(rg$slc45a2$analysis$oDG$data.groups$modelII, add=T,xlim=c(0,3500))

## BNC2

# Do model selection based on the AICc scores
print(rg$bnc2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$bnc2$analysis$oDG));
# Print out the model with the minimum AICc score
print(rg$bnc2$analysis$model.name <-
        rownames(rg$bnc2$analysis$AICcTable
        )[[ which.min(rg$bnc2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rg$bnc2$analysis$model.selected <-
  rg$bnc2$analysis$oDG$data.groups[[rg$bnc2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$bnc2$analysis$model.selected,
                         names(rg$bnc2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$bnc2$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$bnc2$analysis$model.selected, add=T,xlim=c(0,3500))

## ICE1-lncRNA

# Do model selection based on the AICc scores
print(rg$ice1.lncrna$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg$ice1.lncrna$analysis$oDG));
# Print out the model with the minimum AICc score
print(rg$ice1.lncrna$analysis$model.name <-
        rownames(rg$ice1.lncrna$analysis$AICcTable
        )[[ which.min(rg$ice1.lncrna$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
rg$ice1.lncrna$analysis$model.selected <-
  rg$ice1.lncrna$analysis$oDG$data.groups[[rg$ice1.lncrna$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg$ice1.lncrna$analysis$model.selected,
                         names(rg$ice1.lncrna$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg$ice1.lncrna$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(rg$ice1.lncrna$analysis$model.selected, add=T,xlim=c(0,3500))


### Summaries for tytleri-gutturalis clines--------------------------------

## KITLG

# Do model selection based on the AICc scores
print(tg$kitlg$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$kitlg$analysis$oDG));
# Print out the model with the minimum AICc score
print(tg$kitlg$analysis$model.name <-
        rownames(tg$kitlg$analysis$AICcTable
        )[[ which.min(tg$kitlg$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
tg$kitlg$analysis$model.selected <-
  tg$kitlg$analysis$oDG$data.groups[[tg$kitlg$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$kitlg$analysis$model.selected,
                         names(tg$kitlg$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$kitlg$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$kitlg$analysis$model.selected, add =T,xlim=c(0,3500))

## SLC45A2

# Do model selection based on the AICc scores
print(tg$slc45a2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$slc45a2$analysis$oDG));
# Print out the model with the minimum AICc score
print(tg$slc45a2$analysis$model.name <-
        rownames(tg$slc45a2$analysis$AICcTable
        )[[ which.min(tg$slc45a2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
tg$slc45a2$analysis$model.selected <-
  tg$slc45a2$analysis$oDG$data.groups[[tg$slc45a2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$slc45a2$analysis$model.selected,
                         names(tg$slc45a2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$slc45a2$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$slc45a2$analysis$model.selected, add=T,xlim=c(0,3500));

## BNC2

# Do model selection based on the AICc scores
print(tg$bnc2$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg$bnc2$analysis$oDG));
# Print out the model with the minimum AICc score
print(tg$bnc2$analysis$model.name <-
        rownames(tg$bnc2$analysis$AICcTable
        )[[ which.min(tg$bnc2$analysis$AICcTable$AICc )]])
# Extract the hzar.dataGroup object for the selected model
tg$bnc2$analysis$model.selected <-
  tg$bnc2$analysis$oDG$data.groups[[tg$bnc2$analysis$model.name]]
# Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg$bnc2$analysis$model.selected,
                         names(tg$bnc2$analysis$model.selected$data.param)));
# Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg$bnc2$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(tg$bnc2$analysis$model.selected, add=T,xlim=c(0,3500))
