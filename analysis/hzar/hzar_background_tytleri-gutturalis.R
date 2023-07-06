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

data.tg <- read.table(paste0('./input_background_hzar/data.tytleri-gutturalis.',args[1],'.txt'),header=T)
names(data.tg)

## Subset to get relevant data
data.tg <- dplyr::select(data.tg, loc, long, lat, h)
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
if(length(apropos("^tg$",ignore.case=FALSE)) == 0 ||
   !is.list(tg) ) tg <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
tg[[locus]] <- list();
## Space to hold the observed data
tg[[locus]]$obs <- list();
## Space to hold the models to fit
tg[[locus]]$models <- list();
## Space to hold the compiled fit requests
tg[[locus]]$fitRs <- list();
## Space to hold the output data chains
tg[[locus]]$runs <- list();
## Space to hold the analysed data
tg[[locus]]$analysis <- list();

## Assign data for genome background
tg[[locus]]$obs <-
  hzar.doMolecularData1DPops(tg.dists$km,
                             tg.dists$h_mean,
                             tg.dists$h_count);

## Make a helper function to set cline models
tg.loadlocusmodel <- function(scaling,tails,
                             id=paste(scaling,tails,sep="."))
  tg[[locus]]$models[[id]] <<- hzar.makeCline1DFreq(tg[[locus]]$obs, scaling, tails)

tg.loadlocusmodel("fixed","none","modelI");
tg.loadlocusmodel("free" ,"none","modelII");
tg.loadlocusmodel("free" ,"both","modelIII");
tg.loadlocusmodel("free" ,"right","modelIV");
tg.loadlocusmodel("free" ,"left","modelV");
tg.loadlocusmodel("free" ,"mirror","modelVI");

## Check the default settings
print(tg[[locus]]$models)

## Modify all models to focus on the transect region.
## Observations were between 0 and 2835.757 km, with Yekatarinbutg.
## Observations were between 0 and 1942.173 km after removing Yekatarinbutg.

min(tg.dists$km)
max(tg.dists$km)

tg[[locus]]$models <- sapply(tg[[locus]]$models,
                         hzar.model.addBoxReq,
                         -30 , 1973,
                         simplify=FALSE)

## Check the updated settings
print(tg[[locus]]$models)

## Compile each of the models to prepare for fitting
tg[[locus]]$fitRs$init <- sapply(tg[[locus]]$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=tg[[locus]]$obs,
                             verbose=FALSE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
tg[[locus]]$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg[[locus]]$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg[[locus]]$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

tg[[locus]]$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg[[locus]]$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg[[locus]]$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

tg[[locus]]$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg[[locus]]$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg[[locus]]$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg[[locus]]$fitRs$init$modelIV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg[[locus]]$fitRs$init$modelIV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg[[locus]]$fitRs$init$modelIV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg[[locus]]$fitRs$init$modelV$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg[[locus]]$fitRs$init$modelV$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg[[locus]]$fitRs$init$modelV$mcmcParam$seed[[1]] <-
  mainSeed$C 

tg[[locus]]$fitRs$init$modelVI$mcmcParam$chainLength <-
  chainLength;                          #1e5
tg[[locus]]$fitRs$init$modelVI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
tg[[locus]]$fitRs$init$modelVI$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
print(tg[[locus]]$fitRs$init)

## Run the models for an initial chain
tg[[locus]]$runs$init <- list()

tg[[locus]]$runs$init$modelI <-
  hzar.doFit(tg[[locus]]$fitRs$init$modelI)

tg[[locus]]$runs$init$modelII <-
  hzar.doFit(tg[[locus]]$fitRs$init$modelII)

tg[[locus]]$runs$init$modelIII <-
  hzar.doFit(tg[[locus]]$fitRs$init$modelIII)

tg[[locus]]$runs$init$modelIV <-
  hzar.doFit(tg[[locus]]$fitRs$init$modelIV)

tg[[locus]]$runs$init$modelV <-
  hzar.doFit(tg[[locus]]$fitRs$init$modelV)

tg[[locus]]$runs$init$modelVI <-
  hzar.doFit(tg[[locus]]$fitRs$init$modelVI)

## Compile a new set of fit requests using the initial chains 
tg[[locus]]$fitRs$chains <-
  lapply(tg[[locus]]$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
tg[[locus]]$fitRs$chains <-
  hzar.multiFitRequest(tg[[locus]]$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Go ahead and run a chain of 3 runs for every fit request
tg[[locus]]$runs$chains <-  hzar.doChain.multi(tg[[locus]]$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
                lapply(tg[[locus]]$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelII converge?
summary(do.call(mcmc.list,
                lapply(tg[[locus]]$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIII converge?
summary(do.call(mcmc.list,
                lapply(tg[[locus]]$runs$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelIV converge?
summary(do.call(mcmc.list,
                lapply(tg[[locus]]$runs$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelV converge?
summary(do.call(mcmc.list,
                lapply(tg[[locus]]$runs$chains[13:15],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did modelVI converge?
summary(do.call(mcmc.list,
                lapply(tg[[locus]]$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## All three models have three convergent chains, so additional runs are not needed.

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
tg[[locus]]$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(tg[[locus]]$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
tg[[locus]]$analysis$initDGs$modelI <-
  hzar.dataGroup.add(tg[[locus]]$runs$init$modelI)
tg[[locus]]$analysis$initDGs$modelII <-
  hzar.dataGroup.add(tg[[locus]]$runs$init$modelII)
tg[[locus]]$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(tg[[locus]]$runs$init$modelIII)
tg[[locus]]$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(tg[[locus]]$runs$init$modelIV)
tg[[locus]]$analysis$initDGs$modelV <-
  hzar.dataGroup.add(tg[[locus]]$runs$init$modelV)
tg[[locus]]$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(tg[[locus]]$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII, model IV, model V, and model VI).
tg[[locus]]$analysis$oDG <-
  hzar.make.obsDataGroup(tg[[locus]]$analysis$initDGs)
tg[[locus]]$analysis$oDG <-
  hzar.copyModelLabels(tg[[locus]]$analysis$initDGs,
                       tg[[locus]]$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
tg[[locus]]$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(tg[[locus]]$runs$chains,
                                hzar.dataGroup.add),
                         tg[[locus]]$analysis$oDG);

## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII, etc. in the hzar.obsDataGroup object.
print(summary(tg[[locus]]$analysis$oDG$data.groups))

## Do model selection based on the AICc scores
print(tg[[locus]]$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(tg[[locus]]$analysis$oDG));

## Print out the model with the minimum AICc score
print(tg[[locus]]$analysis$model.name <-
        rownames(tg[[locus]]$analysis$AICcTable
        )[[ which.min(tg[[locus]]$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
tg[[locus]]$analysis$model.selected <-
  tg[[locus]]$analysis$oDG$data.groups[[tg[[locus]]$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(tg[[locus]]$analysis$model.selected,
                         names(tg[[locus]]$analysis$model.selected$data.param)));

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(tg[[locus]]$analysis$model.selected))

## Save results to file
save(tg, file = paste0('./results_hzar/hzar_tytleri-gutturalis.background_',args[1],'.RData'))

## End Analysis

quit(save="no")
