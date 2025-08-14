#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Note: this script expects a command line argument in the form of 'locus1', for example.

## Load libraries

#install.packages("hzar") --> had to install from source because HZAR has been removed from CRAN(?)
#install.packages('./hzar_0.2-5.tar.gz', repos = NULL, type ="source")

library(tidyverse)
library(reshape2)
library(hzar)
library(doMC)

## Read in and format data for HZAR (rustica-gutturalis)----------------------

# These data include the sample, population, lat, long, locality, hybrid index
# (background), and hybrid index (outlier regions).

data.rg <- read.table(paste0('./input_background_hzar/data.rustica-gutturalis.',args[1],'.txt'),header=T)
names(data.rg)

## Subset to get relevant data
data.rg <- dplyr::select(data.rg, loc, long, lat, h)
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

## Calculate distances (km) for transect
rg.dists<-hzar.map.latLongSites(data.rg.hzar.cast$loc, data.rg.hzar.cast$lat_mean, data.rg.hzar.cast$long_mean, degrees = TRUE)

join_func<-function(y){
  y<-left_join(y, data.rg.hzar.cast, by=c("site"="loc"))
  y<-y[,-c(2:3)]
  y<-arrange(y, km)
}

## Subset points to the sampling transects
rg.dists<-filter(rg.dists, site=="karasuk" | site=="urumqi" | site=="wuwei" | site=="Lanzhou" | site=="Jiuquan" | site=="Zhangye" | site=="gaotai" | site=="yumen" | site=="xian" | site=="zhengzhou") #urumqi the straight line

##Add km transect distances
rg.dists$km <- hzar.map.distanceFromSite(rg.dists,"karasuk",units="Km")
rg.dists<-join_func(rg.dists)
print(rg.dists)

## Set up MCMC chain parameters------------------------------------------

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

### Analysis on background locus--------------------------------------------

## Set locus name for analysis

locus <- c(args[1])

## Blank out space in memory to hold molecular analysis
if(length(apropos("^rg$",ignore.case=FALSE)) == 0 ||
   !is.list(rg) ) rg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rg[[locus]] <- list();
## Space to hold the observed data
rg[[locus]]$obs <- list();
## Space to hold the models to fit
rg[[locus]]$models <- list();
## Space to hold the compiled fit requests
rg[[locus]]$fitRs <- list();
## Space to hold the output data chains
rg[[locus]]$runs <- list();
## Space to hold the analysed data
rg[[locus]]$analysis <- list();

## Assign data for genome background
rg[[locus]]$obs <-
  hzar.doMolecularData1DPops(rg.dists$km,
                             rg.dists$h_mean,
                             rg.dists$h_count);

## Make a helper function to set cline models
rg.loadlocusmodel <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  rg[[locus]]$models[[id]] <<- hzar.makeCline1DFreq(rg[[locus]]$obs, scaling, tails)

rg.loadlocusmodel("fixed","none","modelI");
rg.loadlocusmodel("free" ,"none","modelII");
rg.loadlocusmodel("free" ,"both","modelIII");
rg.loadlocusmodel("free" ,"right","modelIV");
rg.loadlocusmodel("free" ,"left","modelV");
rg.loadlocusmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rg[[locus]]$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rg.dists$km)
max(rg.dists$km)

rg[[locus]]$models <- sapply(rg[[locus]]$models,
                         hzar.model.addBoxReq,
                         -30 , 1973,
                         simplify=FALSE)

## Check the updated settings
print(rg[[locus]]$models)

## Compile each of the models to prepare for fitting
rg[[locus]]$fitRs$init <- sapply(rg[[locus]]$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=rg[[locus]]$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
rg[[locus]]$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg[[locus]]$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg[[locus]]$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rg[[locus]]$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg[[locus]]$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg[[locus]]$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rg[[locus]]$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg[[locus]]$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg[[locus]]$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg[[locus]]$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg[[locus]]$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg[[locus]]$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg[[locus]]$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg[[locus]]$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg[[locus]]$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rg[[locus]]$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rg[[locus]]$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rg[[locus]]$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rg[[locus]]$fitRs$init)

## Run the models for an initial chain
rg[[locus]]$runs$init <- list()

rg[[locus]]$runs$init$modelI <-
  hzar.doFit(rg[[locus]]$fitRs$init$modelI)

rg[[locus]]$runs$init$modelII <-
  hzar.doFit(rg[[locus]]$fitRs$init$modelII)

rg[[locus]]$runs$init$modelIII <-
  hzar.doFit(rg[[locus]]$fitRs$init$modelIII)

rg[[locus]]$runs$init$modelIV <-
  hzar.doFit(rg[[locus]]$fitRs$init$modelIV)

rg[[locus]]$runs$init$modelV <-
  hzar.doFit(rg[[locus]]$fitRs$init$modelV)

rg[[locus]]$runs$init$modelVI <-
  hzar.doFit(rg[[locus]]$fitRs$init$modelVI)

## Compile a new set of fit requests using the initial chains 
rg[[locus]]$fitRs$chains <-
  lapply(rg[[locus]]$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rg[[locus]]$fitRs$chains <-
  hzar.multiFitRequest(rg[[locus]]$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rg[[locus]]$runs$chains <-  hzar.doChain.multi(rg[[locus]]$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rg[[locus]]$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rg[[locus]]$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rg[[locus]]$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rg[[locus]]$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rg[[locus]]$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rg[[locus]]$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rg[[locus]]$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rg[[locus]]$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rg[[locus]]$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rg[[locus]]$runs$init$modelI)
rg[[locus]]$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rg[[locus]]$runs$init$modelII)
rg[[locus]]$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rg[[locus]]$runs$init$modelIII)
rg[[locus]]$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rg[[locus]]$runs$init$modelIV)
rg[[locus]]$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rg[[locus]]$runs$init$modelV)
rg[[locus]]$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rg[[locus]]$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rg[[locus]]$analysis$oDG <-
  hzar.make.obsDataGroup(rg[[locus]]$analysis$initDGs)
rg[[locus]]$analysis$oDG <-
  hzar.copyModelLabels(rg[[locus]]$analysis$initDGs,
                       rg[[locus]]$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rg[[locus]]$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rg[[locus]]$runs$chains,
                                hzar.dataGroup.add),
                         rg[[locus]]$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rg[[locus]]$analysis$oDG$data.groups))

## Do model selection based on the AICc scores
print(rg[[locus]]$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rg[[locus]]$analysis$oDG));

## Print out the model with the minimum AICc score
print(rg[[locus]]$analysis$model.name <-
        rownames(rg[[locus]]$analysis$AICcTable
        )[[ which.min(rg[[locus]]$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rg[[locus]]$analysis$model.selected <-
  rg[[locus]]$analysis$oDG$data.groups[[rg[[locus]]$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rg[[locus]]$analysis$model.selected,
                         names(rg[[locus]]$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rg[[locus]]$analysis$model.selected))

## Save results to file
save(rg, file = paste0('./results_hzar/hzar_rustica-gutturalis.background_',args[1],'.RData'))

## End Analysis

quit(save="no")
