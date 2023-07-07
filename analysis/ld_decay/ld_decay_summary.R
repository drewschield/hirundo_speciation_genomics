#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Note: this script expects two command line argument of a population name
## in the form of 'rustica', for example, and either 'auto' or 'chrZ'.

### Working directory and dependencies--------------------------------------

library(data.table)
library(tidyverse)

### Read in data-------------------------------------------------------------

r2 <- read.delim(paste0("./summary/pop.",args[1],".",args[2],".hap.ld.summary"),sep="",header=F,check.names=F,stringsAsFactors=F)

### Calculate mean and quantile values with physical distance----------------

colnames(r2) <- c("dist","rsq")

distances <- c()
means <- c()
medians <- c()
q1s <- c()
q3s <- c()
## Loop to calculate summary r2 stats in bins of 100 bp distance, up to 25 kb
test.range <- seq(0,24900,100) # to make range with intervals of 100
for(i in test.range){
  e <- i + 100
  rsqs <- r2 %>% filter(dist >= i & dist < e)
  sums <- summary(rsqs$rsq)
  m <- unname(sums[4]) # unname gets rid of string header in summary output
  d <- unname(sums[3])
  l <- unname(quantile(rsqs$rsq,0.0225))
  u <- unname(quantile(rsqs$rsq,0.975))
  distances=append(distances,e)
  means=append(means,m)
  medians=append(medians,d)
  q1s=append(q1s,l)
  q3s=append(q3s,u)
  print(paste("yeehaw",m))
}

## Format data table and save results to file--------------------------------

decay <- data.frame(distances,means,medians,q1s,q3s)
colnames(decay) <- c("distance","mean","median","Q1","Q3")
write.table(decay,file=paste0("./summary/pop.",args[1],".",args[2],".hap.ld.decay.txt"),row.names=F,quote=F,sep="\t")

## End Analysis

quit(save="no")