############################################################################
# Barn swallow geographic clines in hybrid zone transects
############################################################################

# This script contains commands to plot the results of geographic cline
# analyses in HZAR.

# This workflow involves:
# 1. Loading HZAR Rdata objects
# 2. Plotting of cline models with statistical summaries for center/width

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

### Plot background loci and candidates------------------------------------

par(mfrow=c(1,3))

plot(1,type='n',xlim=c(0,3500),ylim=c(0,1),ylab='Hybrid Index',xlab='Distance (km)')
load(file='./results_hzar/hzar_rustica-tytleri.background_locus1.RData')
hzar.plot.cline(rt[[1]]$analysis$model.selected,pch=NA,col=alpha('grey',0.25),xlim=c(0,2000),add=T)
for (i in rt.list$V1){
  print(i)
  load(file=paste0('./results_hzar/hzar_rustica-tytleri.background_',i,'.RData'))
  hzar.plot.cline(rt[[1]]$analysis$model.selected,add=T,pch=NA,col=alpha('grey',0.25),xlim=c(0,2000))
}

hzar.plot.cline(rt.cand$kitlg$analysis$model.selected,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$plxnc1$analysis$model.selected,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$slc45a2$analysis$oDG$data.groups$modelII,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$bnc2$analysis$model.selected,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$gnaq$analysis$model.selected,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$ror2$analysis$oDG$data.groups$modelII,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$apc.camk4$analysis$oDG$data.groups$modelII,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$ice1.lncrna$analysis$model.selected,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)
hzar.plot.cline(rt.cand$pde1c$analysis$oDG$data.groups$modelII,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)

plot(1,type='n',xlim=c(0,3500),ylim=c(0,1),ylab='Hybrid Index',xlab='Distance (km)')
load(file='./results_hzar/hzar_rustica-gutturalis.background_locus1.RData')
hzar.plot.cline(rg[[1]]$analysis$model.selected,pch=NA,col=alpha('grey',0.25),add=T)
for (i in rg.list$V1){
  print(i)
  load(file=paste0('./results_hzar/hzar_rustica-gutturalis.background_',i,'.RData'))
  hzar.plot.cline(rg[[1]]$analysis$model.selected,add=T,pch=NA,col=alpha('grey',0.25),xlim=c(0,3500))
}

hzar.plot.cline(rg.cand$kitlg$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$plxnc1$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$spef2.prlr$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$slc45a2$analysis$oDG$data.groups$modelII,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$bnc2$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$gnaq$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$ror2$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$apc.camk4$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)
hzar.plot.cline(rg.cand$pde1c$analysis$model.selected,pch=NA,col='purple',add=T,xlim=c(0,3500),lwd=1.5)

plot(1,type='n',xlim=c(0,3500),ylim=c(0,1),ylab='Hybrid Index',xlab='Distance (km)')
load(file='./results_hzar/hzar_tytleri-gutturalis.background_locus1.RData')
hzar.plot.cline(tg[[1]]$analysis$model.selected,pch=NA,col=alpha('grey',0.25),add=T,xlim=c(0,2200))
for (i in tg.list$V1){
  print(i)
  load(file=paste0('./results_hzar/hzar_tytleri-gutturalis.background_',i,'.RData'))
  hzar.plot.cline(tg[[1]]$analysis$model.selected,add=T,pch=NA,col=alpha('grey',0.25),xlim=c(0,2200))
}

hzar.plot.cline(tg.cand$kitlg$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$plxnc1$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$spef2.prlr$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$slc45a2$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$bnc2$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$gnaq$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$ror2$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$apc.camk4$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$ice1.lncrna$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)
hzar.plot.cline(tg.cand$pde1c$analysis$model.selected,pch=NA,col='seagreen',add=T,xlim=c(0,2200),lwd=1.5)

### Read in background cline parameter data---

rt.back <- read.table('./background_locus.param.rustica-tytleri.txt',header=T)
rg.back <- read.table('./background_locus.param.rustica-gutturalis.txt',header=T)
tg.back <- read.table('./background_locus.param.tytleri-gutturalis.txt',header=T)

### Plot background & candidate parameters----

## Cline centers & widths

plot(1,type='n',xlim=c(0,3500),ylim=c(0,12))
segments(x0=mean(rt.back$center)-((mean(rt.back$width))/2),y0=12,x1=mean(rt.back$center)+((mean(rt.back$width))/2),y1=12,col='grey',lwd=1.5)
points(x=mean(rt.back$center),pch=20,y=12,col='grey')
segments(x0=rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$center-((rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y0=11,x1=rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$center+((rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y1=11,col='darkorange',lwd=1.5)
points(x=rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$center,pch=20,y=11,col='darkorange')
segments(x0=rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center-((rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),y0=10,x1=rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center+((rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),y1=10,col='darkorange',lwd=1.5)
points(x=rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center,pch=20,y=10,col='darkorange')
segments(x0=rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=9,x1=rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=9,col='darkorange',lwd=1.5)
points(x=rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=9,col='darkorange')
segments(x0=rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=8,x1=rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=8,col='darkorange',lwd=1.5)
points(x=rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=8,col='darkorange')
segments(x0=rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$center-((rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y0=7,x1=rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$center+((rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y1=7,col='darkorange3',lwd=1.5)
points(x=rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=7,col='darkorange3')
segments(x0=rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$center-((rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)/2),y0=6,x1=rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$center+((rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)/2),y1=6,col='darkorange3',lwd=1.5)
points(x=rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$center,pch=20,y=6,col='darkorange3')
segments(x0=rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=5,x1=rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=5,col='darkorange',lwd=1.5)
points(x=rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=5,col='darkorange')
segments(x0=rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=4,x1=rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=4,col='darkorange',lwd=1.5)
points(x=rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=4,col='darkorange')
segments(x0=rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center-((rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$width)/2),y0=3,x1=rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center+((rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$width)/2),y1=3,col='darkorange',lwd=1.5)
points(x=rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center,pch=20,y=3,col='darkorange')
segments(x0=rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$center-((rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$width)/2),y0=2,x1=rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$center+((rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$width)/2),y1=2,col='darkorange',lwd=1.5)
points(x=rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$center,pch=20,y=2,col='darkorange')
segments(x0=rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=1,x1=rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=1,col='darkorange',lwd=1.5)
points(x=rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=1,col='darkorange')

plot(1,type='n',xlim=c(0,3500),ylim=c(0,12))
segments(x0=mean(rg.back$center)-((mean(rg.back$width))/2),y0=12,x1=mean(rg.back$center)+((mean(rg.back$width))/2),y1=12,col='grey',lwd=1.5)
points(x=mean(rg.back$center),pch=20,y=12,col='grey')
segments(x0=rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center-((rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y0=11,x1=rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center+((rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y1=11,col='purple',lwd=1.5)
points(x=rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center,pch=20,y=11,col='purple')
segments(x0=rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center-((rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),y0=10,x1=rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center+((rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),y1=10,col='purple',lwd=1.5)
points(x=rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center,pch=20,y=10,col='purple')
segments(x0=rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center-((rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)/2),y0=9,x1=rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center+((rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)/2),y1=9,col='purple',lwd=1.5)
points(x=rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center,pch=20,y=9,col='purple')
segments(x0=rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=8,x1=rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=8,col='purple',lwd=1.5)
points(x=rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=8,col='purple')
segments(x0=rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center-((rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y0=7,x1=rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center+((rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y1=7,col='purple3',lwd=1.5)
points(x=rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=7,col='purple3')
segments(x0=rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center-((rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)/2),y0=6,x1=rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center+((rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)/2),y1=6,col='purple3',lwd=1.5)
points(x=rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center,pch=20,y=6,col='purple3')
segments(x0=rg.cand$ror2$analysis$model.selected$ML.cline$param.all$center-((rg.cand$ror2$analysis$model.selected$ML.cline$param.all$width)/2),y0=5,x1=rg.cand$ror2$analysis$model.selected$ML.cline$param.all$center+((rg.cand$ror2$analysis$model.selected$ML.cline$param.all$width)/2),y1=5,col='purple',lwd=1.5)
points(x=rg.cand$ror2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=5,col='purple')
segments(x0=rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center-((rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$width)/2),y0=4,x1=rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center+((rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$width)/2),y1=4,col='purple',lwd=1.5)
points(x=rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center,pch=20,y=4,col='purple')
segments(x0=rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=3,x1=rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=3,col='purple',lwd=1.5)
points(x=rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=3,col='purple')
segments(x0=rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$center-((rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$width)/2),y0=2,x1=rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$center+((rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$width)/2),y1=2,col='purple',lwd=1.5)
points(x=rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$center,pch=20,y=2,col='purple')
segments(x0=rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center-((rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$width)/2),y0=1,x1=rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center+((rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$width)/2),y1=1,col='purple',lwd=1.5)
points(x=rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center,pch=20,y=1,col='purple')

plot(1,type='n',xlim=c(0,3500),ylim=c(0,12))
segments(x0=mean(tg.back$center)-((mean(tg.back$width))/2),y0=12,x1=mean(tg.back$center)+((mean(tg.back$width))/2),y1=12,col='grey',lwd=1.5)
points(x=mean(tg.back$center),pch=20,y=12,col='grey')
segments(x0=tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center-((tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y0=11,x1=tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center+((tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y1=11,col='seagreen',lwd=1.5)
points(x=tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$center,pch=20,y=11,col='seagreen')
segments(x0=tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center-((tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),y0=10,x1=tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center+((tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$width)/2),y1=10,col='seagreen',lwd=1.5)
points(x=tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$center,pch=20,y=10,col='seagreen')
segments(x0=tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center-((tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)/2),y0=9,x1=tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center+((tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$width)/2),y1=9,col='seagreen',lwd=1.5)
points(x=tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$center,pch=20,y=9,col='seagreen')
segments(x0=tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$center-((tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),y0=8,x1=tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$center+((tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),y1=8,col='seagreen',lwd=1.5)
points(x=tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=8,col='seagreen')
segments(x0=tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center-((tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y0=7,x1=tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center+((tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y1=7,col='seagreen3',lwd=1.5)
points(x=tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=7,col='seagreen3')
segments(x0=tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center-((tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)/2),y0=6,x1=tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center+((tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$width)/2),y1=6,col='seagreen3',lwd=1.5)
points(x=tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$center,pch=20,y=6,col='seagreen3')
segments(x0=tg.cand$ror2$analysis$model.selected$ML.cline$param.all$center-((tg.cand$ror2$analysis$model.selected$ML.cline$param.all$width)/2),y0=5,x1=tg.cand$ror2$analysis$model.selected$ML.cline$param.all$center+((tg.cand$ror2$analysis$model.selected$ML.cline$param.all$width)/2),y1=5,col='seagreen',lwd=1.5)
points(x=tg.cand$ror2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=5,col='seagreen')
segments(x0=tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center-((tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$width)/2),y0=4,x1=tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center+((tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$width)/2),y1=4,col='seagreen',lwd=1.5)
points(x=tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$center,pch=20,y=4,col='seagreen')
segments(x0=tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center-((tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$width)/2),y0=3,x1=tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center+((tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$width)/2),y1=3,col='seagreen',lwd=1.5)
points(x=tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$center,pch=20,y=3,col='seagreen')
segments(x0=tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$center-((tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$width)/2),y0=2,x1=tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$center+((tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$width)/2),y1=2,col='seagreen',lwd=1.5)
points(x=tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$center,pch=20,y=2,col='seagreen')
segments(x0=tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center-((tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$width)/2),y0=1,x1=tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center+((tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$width)/2),y1=1,col='seagreen',lwd=1.5)
points(x=tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$center,pch=20,y=1,col='seagreen')

## Output at 3 x 15

## Cline delta h (hybrid index range)
plot(1,type='n',xlim=c(0,12),ylim=c(0,1))
segments(y0=mean(rt.back$pmin),x0=1,y1=mean(rt.back$pmax),x1=1,col='grey',lwd=1.5)

segments(y0=rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMin,x0=2,y1=rt.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMax,x1=2,col='darkorange',lwd=1.5)
segments(y0=rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMin,x0=3,y1=rt.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMax,x1=3,col='darkorange',lwd=1.5)
segments(y0=rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=4,y1=rt.cand$spef2.prlr$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=4,col='darkorange',lwd=1.5)
segments(y0=rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=5,y1=rt.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=5,col='darkorange',lwd=1.5)
segments(y0=rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMin,x0=6,y1=rt.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMax,x1=6,col='darkorange3',lwd=1.5)
segments(y0=rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMin,x0=7,y1=rt.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMax,x1=7,col='darkorange3',lwd=1.5)
segments(y0=rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=8,y1=rt.cand$ror2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=8,col='darkorange',lwd=1.5)
segments(y0=rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=9,y1=rt.cand$apc.camk4$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=9,col='darkorange',lwd=1.5)
segments(y0=rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMin,x0=10,y1=rt.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMax,x1=10,col='darkorange',lwd=1.5)
segments(y0=rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$pMin,x0=11,y1=rt.cand$adcy1$analysis$model.selected$ML.cline$param.all$pMax,x1=11,col='darkorange',lwd=1.5)
segments(y0=rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=12,y1=rt.cand$pde1c$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=12,col='darkorange',lwd=1.5)

plot(1,type='n',xlim=c(0,12),ylim=c(0,1))
segments(y0=mean(rg.back$pmin),x0=1,y1=mean(rg.back$pmax),x1=1,col='grey',lwd=1.5)

segments(y0=rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMin,x0=2,y1=rg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMax,x1=2,col='purple',lwd=1.5)
segments(y0=rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMin,x0=3,y1=rg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMax,x1=3,col='purple',lwd=1.5)
segments(y0=rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMin,x0=4,y1=rg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMax,x1=4,col='purple',lwd=1.5)
segments(y0=rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=5,y1=rg.cand$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=5,col='purple',lwd=1.5)
segments(y0=rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMin,x0=6,y1=rg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMax,x1=6,col='purple3',lwd=1.5)
segments(y0=rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMin,x0=7,y1=rg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMax,x1=7,col='purple3',lwd=1.5)
segments(y0=rg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMin,x0=8,y1=rg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMax,x1=8,col='purple',lwd=1.5)
segments(y0=rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMin,x0=9,y1=rg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMax,x1=9,col='purple',lwd=1.5)
segments(y0=rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMin,x0=10,y1=rg.cand$ice1.lncrna$analysis$oDG$data.groups$modelII$ML.cline$param.all$pMax,x1=10,col='purple',lwd=1.5)
segments(y0=rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$pMin,x0=11,y1=rg.cand$adcy1$analysis$model.selected$ML.cline$param.all$pMax,x1=11,col='purple',lwd=1.5)
segments(y0=rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMin,x0=12,y1=rg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMax,x1=12,col='purple',lwd=1.5)

plot(1,type='n',xlim=c(0,12),ylim=c(0,1))
segments(y0=mean(tg.back$pmin),x0=1,y1=mean(tg.back$pmax),x1=1,col='grey',lwd=1.5)

segments(y0=tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMin,x0=2,y1=tg.cand$kitlg$analysis$model.selected$ML.cline$param.all$pMax,x1=2,col='seagreen',lwd=1.5)
segments(y0=tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMin,x0=3,y1=tg.cand$plxnc1$analysis$model.selected$ML.cline$param.all$pMax,x1=3,col='seagreen',lwd=1.5)
segments(y0=tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMin,x0=4,y1=tg.cand$spef2.prlr$analysis$model.selected$ML.cline$param.all$pMax,x1=4,col='seagreen',lwd=1.5)
segments(y0=tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$pMin,x0=5,y1=tg.cand$slc45a2$analysis$model.selected$ML.cline$param.all$pMax,x1=5,col='seagreen',lwd=1.5)
segments(y0=tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMin,x0=6,y1=tg.cand$bnc2$analysis$model.selected$ML.cline$param.all$pMax,x1=6,col='seagreen3',lwd=1.5)
segments(y0=tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMin,x0=7,y1=tg.cand$gnaq$analysis$model.selected$ML.cline$param.all$pMax,x1=7,col='seagreen3',lwd=1.5)
segments(y0=tg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMin,x0=8,y1=tg.cand$ror2$analysis$model.selected$ML.cline$param.all$pMax,x1=8,col='seagreen',lwd=1.5)
segments(y0=tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMin,x0=9,y1=tg.cand$apc.camk4$analysis$model.selected$ML.cline$param.all$pMax,x1=9,col='seagreen',lwd=1.5)
segments(y0=tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMin,x0=10,y1=tg.cand$ice1.lncrna$analysis$model.selected$ML.cline$param.all$pMax,x1=10,col='seagreen',lwd=1.5)
segments(y0=tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$pMin,x0=11,y1=tg.cand$adcy1$analysis$model.selected$ML.cline$param.all$pMax,x1=11,col='seagreen',lwd=1.5)
segments(y0=tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMin,x0=12,y1=tg.cand$pde1c$analysis$model.selected$ML.cline$param.all$pMax,x1=12,col='seagreen',lwd=1.5)

## Output at 4 x 8


### Plot model I and II clines for loci with similar AIC support----------

par(mfrow=c(4,2))
## rustica-tytleri
hzar.plot.cline(rt$spef2.prlr$analysis$model.selected,col='grey',lty=2,main='rustica-tytleri; SPEF2-PRLR',ylab='h')
hzar.plot.cline(rt$spef2.prlr$analysis$oDG$data.groups$modelII,add=T)

hzar.plot.cline(rt$slc45a2$analysis$model.selected,col='grey',lty=2,main='rustica-tytleri; SLC45A2',ylab='h')
hzar.plot.cline(rt$slc45a2$analysis$oDG$data.groups$modelII,add=T)

hzar.plot.cline(rt$ror2$analysis$model.selected,col='grey',lty=2,main='rustica-tytleri; ROR2',ylab='h')
hzar.plot.cline(rt$ror2$analysis$oDG$data.groups$modelII,add=T)

hzar.plot.cline(rt$apc.camk4$analysis$oDG$data.groups$modelI,col='grey',lty=2,main='rustica-tytleri; APC-CAMK4',ylab='h')
hzar.plot.cline(rt$apc.camk4$analysis$model.selected,add=T)

hzar.plot.cline(rt$pde1c$analysis$model.selected,col='grey',lty=2,main='rustica-tytleri; PDE1C',ylab='h')
hzar.plot.cline(rt$pde1c$analysis$oDG$data.groups$modelII,add=T)

plot(1,type='n',xlim=c(0,2000),ylim=c(0,1))

hzar.plot.cline(rg$slc45a2$analysis$model.selected,col='grey',lty=2,main='rustica-gutturalis; SLC45A2',ylab='h')
hzar.plot.cline(rg$slc45a2$analysis$oDG$data.groups$modelII,add=T)

hzar.plot.cline(rg$ice1.lncrna$analysis$model.selected,col='grey',lty=2,main='rustica-gutturalis; ICE1-lncRNA',ylab='h')
hzar.plot.cline(rg$ice1.lncrna$analysis$oDG$data.groups$modelII,add=T)

### Plot simple clines with dot-whisker plots of center/width--------------

## This compares some top candidates to the genome-wide background estimate

par(mfrow=c(2,3))
## rustica-tytleri
plot(1,type='n',xlim=c(0,3500),ylim=c(0,1),ylab='Hybrid Index',xlab='Distance (km)')
hzar.plot.cline(rt$back$analysis$model.selected,lty=2,lwd=2,xlim=c(0,2000),ylab='Hybrid Index',xlab='Distance (km)',add=T)
hzar.plot.cline(rt$kitlg$analysis$model.selected, add=T,lwd=1.5,xlim=c(0,2000),col='darkorange')
hzar.plot.cline(rt$slc45a2$analysis$oDG$data.groups$modelII, add=T,lwd=1.5,xlim=c(0,2000),col='darkorange2')
hzar.plot.cline(rt$bnc2$analysis$model.selected, add=T,lwd=1.5,xlim=c(0,2000),col='darkorange3')

## rustica-gutturalis
plot(1,type='n',xlim=c(0,3500),ylim=c(0,1),ylab='Hybrid Index',xlab='Distance (km)')
hzar.plot.cline(rg$back$analysis$model.selected,lty=2,lwd=2,xlim=c(0,3500),ylab='Hybrid Index',xlab='Distance (km)',add=T)
hzar.plot.cline(rg$kitlg$analysis$model.selected, add=T,lwd=1.5,xlim=c(0,3500),col='purple')
hzar.plot.cline(rg$slc45a2$analysis$oDG$data.groups$modelII, add=T,lwd=1.5,xlim=c(0,3500),col='purple2')
hzar.plot.cline(rg$bnc2$analysis$model.selected, add=T,lwd=1.5,xlim=c(0,3500),col='purple3')

## tytleri-gutturalis
plot(1,type='n',xlim=c(0,3500),ylim=c(0,1),ylab='Hybrid Index',xlab='Distance (km)')
hzar.plot.cline(tg$back$analysis$model.selected,lty=2,lwd=2,xlim=c(0,2200),ylab='Hybrid Index',xlab='Distance (km)',add=T)
hzar.plot.cline(tg$kitlg$analysis$model.selected, add =T,lwd=1.5,xlim=c(0,2200),col='aquamarine')
hzar.plot.cline(tg$slc45a2$analysis$model.selected, add=T,lwd=1.5,xlim=c(0,2200),col='aquamarine2');
hzar.plot.cline(tg$bnc2$analysis$model.selected, add=T,lwd=1.5,xlim=c(0,2200),col='aquamarine3')

plot(1,type='n',xlim=c(0,3500),ylim=c(0,5))
segments(x0=rt$back$analysis$model.selected$ML.cline$param.all$center-((rt$back$analysis$model.selected$ML.cline$param.all$width)/2),y0=4,x1=rt$back$analysis$model.selected$ML.cline$param.all$center+((rt$back$analysis$model.selected$ML.cline$param.all$width)/2),y1=4,lwd=1.5)
points(x=rt$back$analysis$model.selected$ML.cline$param.all$center,pch=20,y=4)
segments(x0=rt$kitlg$analysis$model.selected$ML.cline$param.all$center-((rt$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y0=3,x1=rt$kitlg$analysis$model.selected$ML.cline$param.all$center+((rt$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y1=3,col='darkorange',lwd=1.5)
points(x=rt$kitlg$analysis$model.selected$ML.cline$param.all$center,pch=20,y=3,col='darkorange')
segments(x0=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=2,x1=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=2,col='darkorange2',lwd=1.5)
points(x=rt$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=2,col='darkorange2')
segments(x0=rt$bnc2$analysis$model.selected$ML.cline$param.all$center-((rt$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y0=1,x1=rt$bnc2$analysis$model.selected$ML.cline$param.all$center+((rt$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y1=1,col='darkorange3',lwd=1.5)
points(x=rt$bnc2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=1,col='darkorange3')

plot(1,type='n',xlim=c(0,3500),ylim=c(0,5))
segments(x0=rg$back$analysis$model.selected$ML.cline$param.all$center-((rg$back$analysis$model.selected$ML.cline$param.all$width)/2),y0=4,x1=rg$back$analysis$model.selected$ML.cline$param.all$center+((rg$back$analysis$model.selected$ML.cline$param.all$width)/2),y1=4,lwd=1.5)
points(x=rg$back$analysis$model.selected$ML.cline$param.all$center,pch=20,y=4)
segments(x0=rg$kitlg$analysis$model.selected$ML.cline$param.all$center-((rg$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y0=3,x1=rg$kitlg$analysis$model.selected$ML.cline$param.all$center+((rg$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y1=3,col='purple',lwd=1.5)
points(x=rg$kitlg$analysis$model.selected$ML.cline$param.all$center,pch=20,y=3,col='purple')
segments(x0=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center-((rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y0=2,x1=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center+((rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$width)/2),y1=2,col='purple2',lwd=1.5)
points(x=rg$slc45a2$analysis$oDG$data.groups$modelII$ML.cline$param.all$center,pch=20,y=2,col='purple2')
segments(x0=rg$bnc2$analysis$model.selected$ML.cline$param.all$center-((rg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y0=1,x1=rg$bnc2$analysis$model.selected$ML.cline$param.all$center+((rg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y1=1,col='purple3',lwd=1.5)
points(x=rg$bnc2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=1,col='purple3')

plot(1,type='n',xlim=c(0,3500),ylim=c(0,5))
segments(x0=tg$back$analysis$model.selected$ML.cline$param.all$center-((tg$back$analysis$model.selected$ML.cline$param.all$width)/2),y0=4,x1=tg$back$analysis$model.selected$ML.cline$param.all$center+((tg$back$analysis$model.selected$ML.cline$param.all$width)/2),y1=4,lwd=1.5)
points(x=tg$back$analysis$model.selected$ML.cline$param.all$center,pch=20,y=4)
segments(x0=tg$kitlg$analysis$model.selected$ML.cline$param.all$center-((tg$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y0=3,x1=tg$kitlg$analysis$model.selected$ML.cline$param.all$center+((tg$kitlg$analysis$model.selected$ML.cline$param.all$width)/2),y1=3,col='aquamarine',lwd=1.5)
points(x=tg$kitlg$analysis$model.selected$ML.cline$param.all$center,pch=20,y=3,col='aquamarine')
segments(x0=tg$slc45a2$analysis$model.selected$ML.cline$param.all$center-((tg$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),y0=2,x1=tg$slc45a2$analysis$model.selected$ML.cline$param.all$center+((tg$slc45a2$analysis$model.selected$ML.cline$param.all$width)/2),y1=2,col='aquamarine2',lwd=1.5)
points(x=tg$slc45a2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=2,col='aquamarine2')
segments(x0=tg$bnc2$analysis$model.selected$ML.cline$param.all$center-((tg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y0=1,x1=tg$bnc2$analysis$model.selected$ML.cline$param.all$center+((tg$bnc2$analysis$model.selected$ML.cline$param.all$width)/2),y1=1,col='aquamarine3',lwd=1.5)
points(x=tg$bnc2$analysis$model.selected$ML.cline$param.all$center,pch=20,y=1,col='aquamarine3')
