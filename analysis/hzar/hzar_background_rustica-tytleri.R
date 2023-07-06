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

## Load data---------------------------------------------------------------
##load("./hzar_background_rustica-tytleri.RData")

## Read in and format data for HZAR (rustica-tytleri)----------------------

# These data include the sample, population, lat, long, locality, hybrid index
# (background), and hybrid index (outlier regions).

data.rt <- read.table(paste0('./input_background_hzar/data.rustica-tytleri.',args[1],'.txt'),header=T)
names(data.rt)

## Subset to get relevant data

data.rt <- dplyr::select(data.rt, loc, long, lat, h)
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
rt.dists$km <- hzar.map.distanceFromSite(rt.dists,"karasuk",units="Km")
rt.dists<-join_func(rt.dists)
rt.dists<-na.omit(rt.dists)
print(rt.dists)

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
if(length(apropos("^rt$",ignore.case=FALSE)) == 0 ||
   !is.list(rt) ) rt <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
rt[[locus]] <- list();
## Space to hold the observed data
rt[[locus]]$obs <- list();
## Space to hold the models to fit
rt[[locus]]$models <- list();
## Space to hold the compiled fit requests
rt[[locus]]$fitRs <- list();
## Space to hold the output data chains
rt[[locus]]$runs <- list();
## Space to hold the analysed data
rt[[locus]]$analysis <- list();

## Assign data for genome background
rt[[locus]]$obs <-
  hzar.doMolecularData1DPops(rt.dists$km,
                             rt.dists$h_mean,
                             rt.dists$h_count);

## Make a helper function to set cline models
rt.loadlocusmodel <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  rt[[locus]]$models[[id]] <<- hzar.makeCline1DFreq(rt[[locus]]$obs, scaling, tails)

rt.loadlocusmodel("fixed","none","modelI");
rt.loadlocusmodel("free" ,"none","modelII");
rt.loadlocusmodel("free" ,"both","modelIII");
rt.loadlocusmodel("free" ,"right","modelIV");
rt.loadlocusmodel("free" ,"left","modelV");
rt.loadlocusmodel("free" ,"mirror","modelVI");

## Check the default settings
print(rt[[locus]]$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinburg.
## Observations were between 0 and 1942.173 km after removing Yekatarinburg.

min(rt.dists$km)
max(rt.dists$km)

rt[[locus]]$models <- sapply(rt[[locus]]$models,
                         hzar.model.addBoxReq,
                         -30 , 1973,
                         simplify=FALSE)

## Check the updated settings
print(rt[[locus]]$models)

## Compile each of the models to prepare for fitting
rt[[locus]]$fitRs$init <- sapply(rt[[locus]]$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=rt[[locus]]$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
rt[[locus]]$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt[[locus]]$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt[[locus]]$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

rt[[locus]]$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt[[locus]]$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt[[locus]]$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

rt[[locus]]$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt[[locus]]$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt[[locus]]$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt[[locus]]$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt[[locus]]$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt[[locus]]$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt[[locus]]$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt[[locus]]$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt[[locus]]$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

rt[[locus]]$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
rt[[locus]]$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
rt[[locus]]$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(rt[[locus]]$fitRs$init)

## Run the models for an initial chain
rt[[locus]]$runs$init <- list()

rt[[locus]]$runs$init$modelI <-
  hzar.doFit(rt[[locus]]$fitRs$init$modelI)

rt[[locus]]$runs$init$modelII <-
  hzar.doFit(rt[[locus]]$fitRs$init$modelII)

rt[[locus]]$runs$init$modelIII <-
  hzar.doFit(rt[[locus]]$fitRs$init$modelIII)

rt[[locus]]$runs$init$modelIV <-
  hzar.doFit(rt[[locus]]$fitRs$init$modelIV)

rt[[locus]]$runs$init$modelV <-
  hzar.doFit(rt[[locus]]$fitRs$init$modelV)

rt[[locus]]$runs$init$modelVI <-
  hzar.doFit(rt[[locus]]$fitRs$init$modelVI)

## Compile a new set of fit requests using the initial chains 
rt[[locus]]$fitRs$chains <-
  lapply(rt[[locus]]$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
rt[[locus]]$fitRs$chains <-
  hzar.multiFitRequest(rt[[locus]]$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
rt[[locus]]$runs$chains <-  hzar.doChain.multi(rt[[locus]]$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(rt[[locus]]$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(rt[[locus]]$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(rt[[locus]]$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(rt[[locus]]$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(rt[[locus]]$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(rt[[locus]]$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
rt[[locus]]$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(rt[[locus]]$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
rt[[locus]]$analysis$initDGs$modelI <-
  hzar.dataGroup.add(rt[[locus]]$runs$init$modelI)
rt[[locus]]$analysis$initDGs$modelII <-
  hzar.dataGroup.add(rt[[locus]]$runs$init$modelII)
rt[[locus]]$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(rt[[locus]]$runs$init$modelIII)
rt[[locus]]$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(rt[[locus]]$runs$init$modelIV)
rt[[locus]]$analysis$initDGs$modelV <-
  hzar.dataGroup.add(rt[[locus]]$runs$init$modelV)
rt[[locus]]$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(rt[[locus]]$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
rt[[locus]]$analysis$oDG <-
  hzar.make.obsDataGroup(rt[[locus]]$analysis$initDGs)
rt[[locus]]$analysis$oDG <-
  hzar.copyModelLabels(rt[[locus]]$analysis$initDGs,
                       rt[[locus]]$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
rt[[locus]]$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(rt[[locus]]$runs$chains,
                                hzar.dataGroup.add),
                         rt[[locus]]$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(rt[[locus]]$analysis$oDG$data.groups))

## Do model selection based on the AICc scores
print(rt[[locus]]$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(rt[[locus]]$analysis$oDG));

## Print out the model with the minimum AICc score
print(rt[[locus]]$analysis$model.name <-
        rownames(rt[[locus]]$analysis$AICcTable
        )[[ which.min(rt[[locus]]$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
rt[[locus]]$analysis$model.selected <-
  rt[[locus]]$analysis$oDG$data.groups[[rt[[locus]]$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(rt[[locus]]$analysis$model.selected,
                         names(rt[[locus]]$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(rt[[locus]]$analysis$model.selected))

## Save results to file
save(rt, file = paste0('./results_hzar/hzar_rustica-tytleri.background_',args[1],'.RData'))

## End Analysis

quit(save="no")
