############################################################################
# Barn swallow genetic diversity, divergence, and differentiation
############################################################################

# We estimated π, dxy, and Fst within and between parental rustica, gutturalis,
# and tytleri populations and their respective hybrid zone using Pixy.

# This script contains commands for parsing, plotting, and performing summary
# statistics on the genome-wide estimates of diversity and differentiation.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/pixy/')

library(data.table)
library(tidyverse)
library(scales)

# options('stringsAsFactors'=FALSE)

### Read in data------------------------------------------------------------

pi.100kb <- read.table('pixy.subspecies.order.pi.100kb.txt',header=T)
pi.10kb <- read.table('pixy.subspecies.order.pi.10kb.txt',header=T)
dxy.100kb <- read.table('pixy.subspecies.order.dxy.100kb.txt',header=T)
dxy.10kb <- read.table('pixy.subspecies.order.dxy.10kb.txt',header=T)
fst.100kb <- read.table('pixy.subspecies.order.fst.100kb.txt',header=T)
fst.10kb <- read.table('pixy.subspecies.order.fst.10kb.txt',header=T)

### Parse populations-------------------------------------------------------

pi.100kb.ru <- pi.100kb %>% filter(str_detect(pop, 'RU'))
pi.100kb.gu <- pi.100kb %>% filter(str_detect(pop, 'GU'))
pi.100kb.ty <- pi.100kb %>% filter(str_detect(pop, 'TY'))
dxy.100kb.ruty <- dxy.100kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.100kb.rugu <- dxy.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'RU'))
dxy.100kb.guty <- dxy.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.100kb.ruty <- fst.100kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.100kb.rugu <- fst.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'RU'))
fst.100kb.guty <- fst.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))

pi.10kb.ru <- pi.10kb %>% filter(str_detect(pop, 'RU'))
pi.10kb.gu <- pi.10kb %>% filter(str_detect(pop, 'GU'))
pi.10kb.ty <- pi.10kb %>% filter(str_detect(pop, 'TY'))
dxy.10kb.ruty <- dxy.10kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
dxy.10kb.rugu <- dxy.10kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'RU'))
dxy.10kb.guty <- dxy.10kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
fst.10kb.ruty <- fst.10kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
fst.10kb.rugu <- fst.10kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'RU'))
fst.10kb.guty <- fst.10kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))

### Convert negative Fst to zero--------------------------------------------

fst.100kb.ruty$avg_wc_fst[fst.100kb.ruty$avg_wc_fst<0] <- 0
fst.100kb.rugu$avg_wc_fst[fst.100kb.rugu$avg_wc_fst<0] <- 0
fst.100kb.guty$avg_wc_fst[fst.100kb.guty$avg_wc_fst<0] <- 0

fst.10kb.ruty$avg_wc_fst[fst.10kb.ruty$avg_wc_fst<0] <- 0
fst.10kb.rugu$avg_wc_fst[fst.10kb.rugu$avg_wc_fst<0] <- 0
fst.10kb.guty$avg_wc_fst[fst.10kb.guty$avg_wc_fst<0] <- 0

### Parse chromosomes-------------------------------------------------------

# W chromosome
pi.100kb.ru.w <- pi.100kb.ru %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
pi.100kb.gu.w <- pi.100kb.gu %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
pi.100kb.ty.w <- pi.100kb.ty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
dxy.100kb.rugu.w <- dxy.100kb.rugu %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
dxy.100kb.ruty.w <- dxy.100kb.ruty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
dxy.100kb.guty.w <- dxy.100kb.guty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
fst.100kb.rugu.w <- fst.100kb.rugu %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
fst.100kb.ruty.w <- fst.100kb.ruty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
fst.100kb.guty.w <- fst.100kb.guty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))

pi.10kb.ru.w <- pi.10kb.ru %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
pi.10kb.gu.w <- pi.10kb.gu %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
pi.10kb.ty.w <- pi.10kb.ty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
dxy.10kb.rugu.w <- dxy.10kb.rugu %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
dxy.10kb.ruty.w <- dxy.10kb.ruty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
dxy.10kb.guty.w <- dxy.10kb.guty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
fst.10kb.rugu.w <- fst.10kb.rugu %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
fst.10kb.ruty.w <- fst.10kb.ruty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))
fst.10kb.guty.w <- fst.10kb.guty %>% filter(str_detect(chromosome,'NC_053487.1') | str_detect(chromosome,'NW_024403836.1') | str_detect(chromosome,'NW_024403837.1'))

# Z chromosome
pi.100kb.ru.z <- pi.100kb.ru %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
pi.100kb.gu.z <- pi.100kb.gu %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
pi.100kb.ty.z <- pi.100kb.ty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
dxy.100kb.rugu.z <- dxy.100kb.rugu %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
dxy.100kb.ruty.z <- dxy.100kb.ruty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
dxy.100kb.guty.z <- dxy.100kb.guty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
fst.100kb.rugu.z <- fst.100kb.rugu %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
fst.100kb.ruty.z <- fst.100kb.ruty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
fst.100kb.guty.z <- fst.100kb.guty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))

pi.10kb.ru.z <- pi.10kb.ru %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
pi.10kb.gu.z <- pi.10kb.gu %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
pi.10kb.ty.z <- pi.10kb.ty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
dxy.10kb.rugu.z <- dxy.10kb.rugu %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
dxy.10kb.ruty.z <- dxy.10kb.ruty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
dxy.10kb.guty.z <- dxy.10kb.guty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
fst.10kb.rugu.z <- fst.10kb.rugu %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
fst.10kb.ruty.z <- fst.10kb.ruty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))
fst.10kb.guty.z <- fst.10kb.guty %>% filter(str_detect(chromosome,'NC_053488.1') | str_detect(chromosome,'NW_024403838.1'))

# Autosomes
pi.100kb.ru.a <- pi.100kb.ru %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
pi.100kb.gu.a <- pi.100kb.gu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
pi.100kb.ty.a <- pi.100kb.ty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
dxy.100kb.rugu.a <- dxy.100kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
dxy.100kb.ruty.a <- dxy.100kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
dxy.100kb.guty.a <- dxy.100kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
fst.100kb.rugu.a <- fst.100kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
fst.100kb.ruty.a <- fst.100kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
fst.100kb.guty.a <- fst.100kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))

pi.10kb.ru.a <- pi.10kb.ru %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
pi.10kb.gu.a <- pi.10kb.gu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
pi.10kb.ty.a <- pi.10kb.ty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
dxy.10kb.rugu.a <- dxy.10kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
dxy.10kb.ruty.a <- dxy.10kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
dxy.10kb.guty.a <- dxy.10kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
fst.10kb.rugu.a <- fst.10kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
fst.10kb.ruty.a <- fst.10kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))
fst.10kb.guty.a <- fst.10kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1') & !str_detect(chromosome,'NC_053488.1') & !str_detect(chromosome,'NW_024403838.1'))

# Autosomes and Z chromosome
pi.100kb.ru <- pi.100kb.ru %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
pi.100kb.gu <- pi.100kb.gu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
pi.100kb.ty <- pi.100kb.ty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
dxy.100kb.rugu <- dxy.100kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
dxy.100kb.ruty <- dxy.100kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
dxy.100kb.guty <- dxy.100kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
fst.100kb.rugu <- fst.100kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
fst.100kb.ruty <- fst.100kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
fst.100kb.guty <- fst.100kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))

pi.10kb.ru <- pi.10kb.ru %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
pi.10kb.gu <- pi.10kb.gu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
pi.10kb.ty <- pi.10kb.ty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
dxy.10kb.rugu <- dxy.10kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
dxy.10kb.ruty <- dxy.10kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
dxy.10kb.guty <- dxy.10kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
fst.10kb.rugu <- fst.10kb.rugu %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
fst.10kb.ruty <- fst.10kb.ruty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))
fst.10kb.guty <- fst.10kb.guty %>% filter(!str_detect(chromosome,'NC_053487.1') & !str_detect(chromosome,'NW_024403836.1') & !str_detect(chromosome,'NW_024403837.1'))

### Calculate population branch statistics (PBS) per subspecies-------------

# Get rid of W chromosome
fst.100kb <- fst.100kb %>% filter(!str_detect(chromosome,'NC_053487.1'))

# Parse Fst results:
er.gu <- fst.100kb %>% filter(str_detect(pop1, 'ER') & str_detect(pop2, 'GU'))
er.ru <- fst.100kb %>% filter(str_detect(pop1, 'ER') & str_detect(pop2, 'RU'))
er.tr <- fst.100kb %>% filter(str_detect(pop1, 'ER') & str_detect(pop2, 'TR'))
er.ty <- fst.100kb %>% filter(str_detect(pop1, 'ER') & str_detect(pop2, 'TY'))
er.sa <- fst.100kb %>% filter(str_detect(pop1, 'ER') & str_detect(pop2, 'SA'))
gu.sa <- fst.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'SA'))
gu.ty <- fst.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'TY'))
ru.gu <- fst.100kb %>% filter(str_detect(pop1, 'GU') & str_detect(pop2, 'RU'))
ru.sa <- fst.100kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'SA'))
ru.tr <- fst.100kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TR'))
ru.ty <- fst.100kb %>% filter(str_detect(pop1, 'RU') & str_detect(pop2, 'TY'))
sa.ty <- fst.100kb %>% filter(str_detect(pop1, 'SA') & str_detect(pop2, 'TY'))

# Merge sets of Fst results and perform analysis:

pbs.er.data <-merge(er.ty,er.sa,by=c("chromosome","window_pos_1"),sort=F) 
pbs.er.data <-merge(pbs.er.data,sa.ty,by=c("chromosome","window_pos_1"),sort=F)
pbs.er <- ((pbs.er.data$avg_wc_fst.x + pbs.er.data$avg_wc_fst.y)-pbs.er.data$avg_wc_fst)/2
pbs.er <-replace(pbs.er, pbs.er<0, 0)
pbs.er.data$pbs.er <- pbs.er

pbs.gu.data <-merge(gu.ty,gu.sa,by=c("chromosome","window_pos_1"),sort=F) 
pbs.gu.data <-merge(pbs.gu.data,sa.ty,by=c("chromosome","window_pos_1"),sort=F)
pbs.gu <- ((pbs.gu.data$avg_wc_fst.x + pbs.gu.data$avg_wc_fst.y)-pbs.gu.data$avg_wc_fst)/2
pbs.gu <-replace(pbs.gu, pbs.gu<0, 0)
pbs.gu.data$pbs.gu <- pbs.gu

pbs.ru.data <-merge(ru.gu,er.ru,by=c("chromosome","window_pos_1"),sort=F) 
pbs.ru.data <-merge(pbs.ru.data,er.gu,by=c("chromosome","window_pos_1"),sort=F)
pbs.ru <- ((pbs.ru.data$avg_wc_fst.x + pbs.ru.data$avg_wc_fst.y)-pbs.ru.data$avg_wc_fst)/2
pbs.ru <-replace(pbs.ru, pbs.ru<0, 0)
pbs.ru.data$pbs.ru <- pbs.ru

pbs.sa.data <-merge(ru.sa,er.sa,by=c("chromosome","window_pos_1"),sort=F) 
pbs.sa.data <-merge(pbs.sa.data,er.ru,by=c("chromosome","window_pos_1"),sort=F)
pbs.sa <- ((pbs.sa.data$avg_wc_fst.x + pbs.sa.data$avg_wc_fst.y)-pbs.sa.data$avg_wc_fst)/2
pbs.sa <-replace(pbs.sa, pbs.sa<0, 0)
pbs.sa.data$pbs.sa <- pbs.sa

pbs.tr.data <-merge(ru.tr,er.tr,by=c("chromosome","window_pos_1"),sort=F) 
pbs.tr.data <-merge(pbs.tr.data,er.ru,by=c("chromosome","window_pos_1"),sort=F)
pbs.tr <- ((pbs.tr.data$avg_wc_fst.x + pbs.tr.data$avg_wc_fst.y)-pbs.tr.data$avg_wc_fst)/2
pbs.tr <-replace(pbs.tr, pbs.tr<0, 0)
pbs.tr.data$pbs.tr <- pbs.tr

pbs.ty.data <-merge(ru.ty,gu.ty,by=c("chromosome","window_pos_1"),sort=F) 
pbs.ty.data <-merge(pbs.ty.data,ru.gu,by=c("chromosome","window_pos_1"),sort=F)
pbs.ty <- ((pbs.ty.data$avg_wc_fst.x + pbs.ty.data$avg_wc_fst.y)-pbs.ty.data$avg_wc_fst)/2
pbs.ty <-replace(pbs.ty, pbs.ty<0, 0)
pbs.ty.data$pbs.ty <- pbs.ty

# Take a look
par(mfrow=c(6,1))
plot(pbs.er.data$pbs.er,type='l')
plot(pbs.gu.data$pbs.gu,type='l')
plot(pbs.ru.data$pbs.ru,type='l')
plot(pbs.sa.data$pbs.sa,type='l')
plot(pbs.tr.data$pbs.tr,type='l')
plot(pbs.ty.data$pbs.ty,type='l')

# Write tables!
write.table(pbs.er.data,file='./pbs/pbs.er.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.gu.data,file='./pbs/pbs.gu.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.ru.data,file='./pbs/pbs.ru.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.sa.data,file='./pbs/pbs.sa.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.tr.data,file='./pbs/pbs.tr.data.txt',quote=F,row.names=F,sep = "\t")
write.table(pbs.ty.data,file='./pbs/pbs.ty.data.txt',quote=F,row.names=F,sep = "\t")


### Calculate statistics in 1 Mb windows-----------------------------------

n <- 10 # every 10 rows

# Autosomes + Z chromosome
setDT(pi.100kb.ru)
pi.1mb.ru <- aggregate(pi.100kb.ru, list(rep(1:(nrow(pi.100kb.ru) %/% n + 1), each = n, len = nrow(pi.100kb.ru))), mean)[-1];
setDT(pi.100kb.gu)
pi.1mb.gu <- aggregate(pi.100kb.gu, list(rep(1:(nrow(pi.100kb.gu) %/% n + 1), each = n, len = nrow(pi.100kb.gu))), mean)[-1];
setDT(pi.100kb.ty)
pi.1mb.ty <- aggregate(pi.100kb.ty, list(rep(1:(nrow(pi.100kb.ty) %/% n + 1), each = n, len = nrow(pi.100kb.ty))), mean)[-1];

setDT(dxy.100kb.rugu)
dxy.1mb.rugu <- aggregate(dxy.100kb.rugu, list(rep(1:(nrow(dxy.100kb.rugu) %/% n + 1), each = n, len = nrow(dxy.100kb.rugu))), mean)[-1];
setDT(dxy.100kb.ruty)
dxy.1mb.ruty <- aggregate(dxy.100kb.ruty, list(rep(1:(nrow(dxy.100kb.ruty) %/% n + 1), each = n, len = nrow(dxy.100kb.ruty))), mean)[-1];
setDT(dxy.100kb.guty)
dxy.1mb.guty <- aggregate(dxy.100kb.guty, list(rep(1:(nrow(dxy.100kb.guty) %/% n + 1), each = n, len = nrow(dxy.100kb.guty))), mean)[-1];

setDT(fst.100kb.rugu)
fst.1mb.rugu <- aggregate(fst.100kb.rugu, list(rep(1:(nrow(fst.100kb.rugu) %/% n + 1), each = n, len = nrow(fst.100kb.rugu))), mean)[-1];
setDT(fst.100kb.ruty)
fst.1mb.ruty <- aggregate(fst.100kb.ruty, list(rep(1:(nrow(fst.100kb.ruty) %/% n + 1), each = n, len = nrow(fst.100kb.ruty))), mean)[-1];
setDT(fst.100kb.guty)
fst.1mb.guty <- aggregate(fst.100kb.guty, list(rep(1:(nrow(fst.100kb.guty) %/% n + 1), each = n, len = nrow(fst.100kb.guty))), mean)[-1];

setDT(pbs.er.data)
pbs.1mb.er <- aggregate(pbs.er.data, list(rep(1:(nrow(pbs.er.data) %/% n + 1), each = n, len = nrow(pbs.er.data))), mean)[-1];
setDT(pbs.gu.data)
pbs.1mb.gu <- aggregate(pbs.gu.data, list(rep(1:(nrow(pbs.gu.data) %/% n + 1), each = n, len = nrow(pbs.gu.data))), mean)[-1];
setDT(pbs.ru.data)
pbs.1mb.ru <- aggregate(pbs.ru.data, list(rep(1:(nrow(pbs.ru.data) %/% n + 1), each = n, len = nrow(pbs.ru.data))), mean)[-1];
setDT(pbs.sa.data)
pbs.1mb.sa <- aggregate(pbs.sa.data, list(rep(1:(nrow(pbs.sa.data) %/% n + 1), each = n, len = nrow(pbs.sa.data))), mean)[-1];
setDT(pbs.tr.data)
pbs.1mb.tr <- aggregate(pbs.tr.data, list(rep(1:(nrow(pbs.tr.data) %/% n + 1), each = n, len = nrow(pbs.tr.data))), mean)[-1];
setDT(pbs.ty.data)
pbs.1mb.ty <- aggregate(pbs.ty.data, list(rep(1:(nrow(pbs.ty.data) %/% n + 1), each = n, len = nrow(pbs.ty.data))), mean)[-1];

# Autosomes
setDT(pi.100kb.ru.a)
pi.1mb.ru.a <- aggregate(pi.100kb.ru.a, list(rep(1:(nrow(pi.100kb.ru.a) %/% n + 1), each = n, len = nrow(pi.100kb.ru.a))), mean)[-1];
setDT(pi.100kb.gu.a)
pi.1mb.gu.a <- aggregate(pi.100kb.gu.a, list(rep(1:(nrow(pi.100kb.gu.a) %/% n + 1), each = n, len = nrow(pi.100kb.gu.a))), mean)[-1];
setDT(pi.100kb.ty.a)
pi.1mb.ty.a <- aggregate(pi.100kb.ty.a, list(rep(1:(nrow(pi.100kb.ty.a) %/% n + 1), each = n, len = nrow(pi.100kb.ty.a))), mean)[-1];

setDT(dxy.100kb.rugu.a)
dxy.1mb.rugu.a <- aggregate(dxy.100kb.rugu.a, list(rep(1:(nrow(dxy.100kb.rugu.a) %/% n + 1), each = n, len = nrow(dxy.100kb.rugu.a))), mean)[-1];
setDT(dxy.100kb.ruty.a)
dxy.1mb.ruty.a <- aggregate(dxy.100kb.ruty.a, list(rep(1:(nrow(dxy.100kb.ruty.a) %/% n + 1), each = n, len = nrow(dxy.100kb.ruty.a))), mean)[-1];
setDT(dxy.100kb.guty.a)
dxy.1mb.guty.a <- aggregate(dxy.100kb.guty.a, list(rep(1:(nrow(dxy.100kb.guty.a) %/% n + 1), each = n, len = nrow(dxy.100kb.guty.a))), mean)[-1];

setDT(fst.100kb.rugu.a)
fst.1mb.rugu.a <- aggregate(fst.100kb.rugu.a, list(rep(1:(nrow(fst.100kb.rugu.a) %/% n + 1), each = n, len = nrow(fst.100kb.rugu.a))), mean)[-1];
setDT(fst.100kb.ruty.a)
fst.1mb.ruty.a <- aggregate(fst.100kb.ruty.a, list(rep(1:(nrow(fst.100kb.ruty.a) %/% n + 1), each = n, len = nrow(fst.100kb.ruty.a))), mean)[-1];
setDT(fst.100kb.guty.a)
fst.1mb.guty.a <- aggregate(fst.100kb.guty.a, list(rep(1:(nrow(fst.100kb.guty.a) %/% n + 1), each = n, len = nrow(fst.100kb.guty.a))), mean)[-1];

# Z chromosome
setDT(pi.100kb.ru.z)
pi.1mb.ru.z <- aggregate(pi.100kb.ru.z, list(rep(1:(nrow(pi.100kb.ru.z) %/% n + 1), each = n, len = nrow(pi.100kb.ru.z))), mean)[-1];
setDT(pi.100kb.gu.z)
pi.1mb.gu.z <- aggregate(pi.100kb.gu.z, list(rep(1:(nrow(pi.100kb.gu.z) %/% n + 1), each = n, len = nrow(pi.100kb.gu.z))), mean)[-1];
setDT(pi.100kb.ty.z)
pi.1mb.ty.z <- aggregate(pi.100kb.ty.z, list(rep(1:(nrow(pi.100kb.ty.z) %/% n + 1), each = n, len = nrow(pi.100kb.ty.z))), mean)[-1];

setDT(dxy.100kb.rugu.z)
dxy.1mb.rugu.z <- aggregate(dxy.100kb.rugu.z, list(rep(1:(nrow(dxy.100kb.rugu.z) %/% n + 1), each = n, len = nrow(dxy.100kb.rugu.z))), mean)[-1];
setDT(dxy.100kb.ruty.z)
dxy.1mb.ruty.z <- aggregate(dxy.100kb.ruty.z, list(rep(1:(nrow(dxy.100kb.ruty.z) %/% n + 1), each = n, len = nrow(dxy.100kb.ruty.z))), mean)[-1];
setDT(dxy.100kb.guty.z)
dxy.1mb.guty.z <- aggregate(dxy.100kb.guty.z, list(rep(1:(nrow(dxy.100kb.guty.z) %/% n + 1), each = n, len = nrow(dxy.100kb.guty.z))), mean)[-1];

setDT(fst.100kb.rugu.z)
fst.1mb.rugu.z <- aggregate(fst.100kb.rugu.z, list(rep(1:(nrow(fst.100kb.rugu.z) %/% n + 1), each = n, len = nrow(fst.100kb.rugu.z))), mean)[-1];
setDT(fst.100kb.ruty.z)
fst.1mb.ruty.z <- aggregate(fst.100kb.ruty.z, list(rep(1:(nrow(fst.100kb.ruty.z) %/% n + 1), each = n, len = nrow(fst.100kb.ruty.z))), mean)[-1];
setDT(fst.100kb.guty.z)
fst.1mb.guty.z <- aggregate(fst.100kb.guty.z, list(rep(1:(nrow(fst.100kb.guty.z) %/% n + 1), each = n, len = nrow(fst.100kb.guty.z))), mean)[-1];


### Recombination rates-----------------------------------------------------

rho <- read.table('../pyrho/rustica.rmap.100kb.txt',header=T)
rho$start.fix <- rho$start+1

### Set chromosome colors---------------------------------------------------

color = rep(NA, length=length(pi.100kb.ru$chromosome))
color[which(pi.100kb.ru$chromosome=="NC_053451.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053453.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053450.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053452.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053454.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053470.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053455.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053457.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053456.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053458.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053459.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053462.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053460.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053461.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053463.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053464.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053466.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053469.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053467.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053468.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053465.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053471.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053477.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053474.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053472.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053478.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053473.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053476.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053475.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053479.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053480.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053481.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053482.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053483.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053484.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053485.1")] = "goldenrod3"
color[which(pi.100kb.ru$chromosome=="NC_053486.1")] = "grey"
color[which(pi.100kb.ru$chromosome=="NC_053488.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053487.1")] = "grey"

#color[which(pi.100kb.ru$chromosome=="NC_053451.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403813.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053453.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403819.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053450.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403809.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403810.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403811.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403812.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053452.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403814.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403815.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403816.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403817.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403818.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053454.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403820.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403821.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403822.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403823.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053470.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403828.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053455.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403824.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053457.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053456.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053458.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053459.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053462.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403826.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053460.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403825.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053461.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053463.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053464.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403827.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053466.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053469.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053467.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053468.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053465.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053471.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403829.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053477.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053474.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053472.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053478.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403832.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403833.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053473.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053476.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403830.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403831.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053475.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053479.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403834.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053480.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053481.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053482.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403835.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053483.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053484.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053485.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053486.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NC_053488.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NW_024403838.1")] = "goldenrod3"
#color[which(pi.100kb.ru$chromosome=="NC_053487.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403836.1")] = "grey"
#color[which(pi.100kb.ru$chromosome=="NW_024403837.1")] = "grey"

### Plot scans on autosomes-------------------------------------------------

par(mfrow=c(3,1))
plot(pi.100kb.ru.a$avg_pi,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='π',main='rustica')
lines(as.integer(row.names(pi.1mb.ru.a))*10,pi.1mb.ru.a$avg_pi,col='black',lwd=1.25)
plot(pi.100kb.ty.a$avg_pi,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='π',main='tytleri')
lines(as.integer(row.names(pi.1mb.ty.a))*10,pi.1mb.ty.a$avg_pi,col='black',lwd=1.25)
plot(pi.100kb.gu.a$avg_pi,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='π',main='gutturalis')
lines(as.integer(row.names(pi.1mb.gu.a))*10,pi.1mb.gu.a$avg_pi,col='black',lwd=1.25)

plot(dxy.100kb.ruty.a$avg_dxy,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='dxy',main='rustica-tytleri')
lines(as.integer(row.names(dxy.1mb.ruty.a))*10,dxy.1mb.ruty.a$avg_dxy,col='black',lwd=1.25)
plot(dxy.100kb.rugu.a$avg_dxy,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='dxy',main='rustica-gutturalis')
lines(as.integer(row.names(dxy.1mb.rugu.a))*10,dxy.1mb.rugu.a$avg_dxy,col='black',lwd=1.25)
plot(dxy.100kb.guty.a$avg_dxy,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='dxy',main='tytleri-gutturalis')
lines(as.integer(row.names(dxy.1mb.guty.a))*10,dxy.1mb.guty.a$avg_dxy,col='black',lwd=1.25)

plot(fst.100kb.ruty.a$avg_wc_fst,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='Fst',main='rustica-tytleri',ylim=c(0,0.4))
lines(as.integer(row.names(fst.1mb.ruty.a))*10,fst.1mb.ruty.a$avg_wc_fst,col='black',lwd=1.25)
plot(fst.100kb.rugu.a$avg_wc_fst,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='Fst',main='rustica-gutturalis',ylim=c(0,0.4))
lines(as.integer(row.names(fst.1mb.rugu.a))*10,fst.1mb.rugu.a$avg_wc_fst,col='black',lwd=1.25)
plot(fst.100kb.guty.a$avg_wc_fst,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='Fst',main='tytleri-gutturalis',ylim=c(0,0.4))
lines(as.integer(row.names(fst.1mb.guty.a))*10,fst.1mb.guty.a$avg_wc_fst,col='black',lwd=1.25)

par(mfrow=c(6,1))
plot(pbs.er.data$pbs.er,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='PBS',main='erythrogaster')
lines(as.integer(row.names(pbs.1mb.er))*10,pbs.1mb.er$pbs.er,col='black',lwd=1.25)
plot(pbs.ty.data$pbs.ty,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='PBS',main='tytleri')
lines(as.integer(row.names(pbs.1mb.ty))*10,pbs.1mb.ty$pbs.ty,col='black',lwd=1.25)
plot(pbs.gu.data$pbs.gu,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='PBS',main='gutturalis')
lines(as.integer(row.names(pbs.1mb.gu))*10,pbs.1mb.gu$pbs.gu,col='black',lwd=1.25)
plot(pbs.ru.data$pbs.ru,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='PBS',main='rustica')
lines(as.integer(row.names(pbs.1mb.ru))*10,pbs.1mb.ru$pbs.ru,col='black',lwd=1.25)
plot(pbs.tr.data$pbs.tr,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='PBS',main='transitiva')
lines(as.integer(row.names(pbs.1mb.tr))*10,pbs.1mb.tr$pbs.tr,col='black',lwd=1.25)
plot(pbs.sa.data$pbs.sa,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='PBS',main='savignii')
lines(as.integer(row.names(pbs.1mb.sa))*10,pbs.1mb.sa$pbs.sa,col='black',lwd=1.25)

### Plot scans on Z chromosome----------------------------------------------

par(mfrow=c(3,1))
plot(pi.100kb.ru.z$avg_pi,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='π',main='rustica')
lines(as.integer(row.names(pi.1mb.ru.z))*10,pi.1mb.ru.z$avg_pi,col='black',lwd=1.25)
plot(pi.100kb.ty.z$avg_pi,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='π',main='tytleri')
lines(as.integer(row.names(pi.1mb.ty.z))*10,pi.1mb.ty.z$avg_pi,col='black',lwd=1.25)
plot(pi.100kb.gu.z$avg_pi,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='π',main='gutturalis')
lines(as.integer(row.names(pi.1mb.gu.z))*10,pi.1mb.gu.z$avg_pi,col='black',lwd=1.25)

plot(dxy.100kb.ruty.z$avg_dxy,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='dxy',main='rustica-tytleri')
lines(as.integer(row.names(dxy.1mb.ruty.z))*10,dxy.1mb.ruty.z$avg_dxy,col='black',lwd=1.25)
plot(dxy.100kb.rugu.z$avg_dxy,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='dxy',main='rustica-gutturalis')
lines(as.integer(row.names(dxy.1mb.rugu.z))*10,dxy.1mb.rugu.z$avg_dxy,col='black',lwd=1.25)
plot(dxy.100kb.guty.z$avg_dxy,pch=20,col=alpha(color,0.1),ylim=c(0,0.015),xlab='Genomic Position (Mb)',ylab='dxy',main='tytleri-gutturalis')
lines(as.integer(row.names(dxy.1mb.guty.z))*10,dxy.1mb.guty.z$avg_dxy,col='black',lwd=1.25)

plot(fst.100kb.ruty.z$avg_wc_fst,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='Fst',main='rustica-tytleri',ylim=c(0,0.8))
lines(as.integer(row.names(fst.1mb.ruty.z))*10,fst.1mb.ruty.z$avg_wc_fst,col='black',lwd=1.25)
plot(fst.100kb.rugu.z$avg_wc_fst,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='Fst',main='rustica-gutturalis',ylim=c(0,0.8))
lines(as.integer(row.names(fst.1mb.rugu.z))*10,fst.1mb.rugu.z$avg_wc_fst,col='black',lwd=1.25)
plot(fst.100kb.guty.z$avg_wc_fst,pch=20,col=alpha(color,0.1),xlab='Genomic Position (Mb)',ylab='Fst',main='tytleri-gutturalis',ylim=c(0,0.8))
lines(as.integer(row.names(fst.1mb.guty.z))*10,fst.1mb.guty.z$avg_wc_fst,col='black',lwd=1.25)

### Plot scans on the W chromosome------------------------------------------

# Plot - output at 6 x 4
par(mfrow=c(3,1))

plot(pi.100kb.ru.w$avg_pi,type='l',col='red',lwd=1.5,ylab='pi',xlab='W chromosome position')
lines(pi.100kb.gu.w$avg_pi,type='l',col='blue',lwd=1.5)
lines(pi.100kb.ty.w$avg_pi,type='l',col='goldenrod',lwd=1.5)

plot(dxy.100kb.ruty.w$avg_dxy,type='l',col='purple',lwd=1.5,ylab='dxy',xlab='W chromosome position')
lines(dxy.100kb.rugu.w$avg_dxy,type='l',col='darkorange',lwd=1.5)
lines(dxy.100kb.guty.w$avg_dxy,type='l',col='seagreen',lwd=1.5)

plot(fst.100kb.ruty.w$avg_wc_fst,type='l',col='purple',lwd=1.5,ylab='Fst',xlab='W chromosome position')
lines(fst.100kb.rugu.w$avg_wc_fst,type='l',col='darkorange',lwd=1.5)
lines(fst.100kb.guty.w$avg_wc_fst,type='l',col='seagreen',lwd=1.5)

### Plot distributions-------------------------------------------------------

par(mfrow=c(1,1))
boxplot(pi.10kb.ru.a$avg_pi,pi.10kb.ty.a$avg_pi,pi.10kb.gu.a$avg_pi,pi.10kb.ru.z$avg_pi,pi.10kb.ty.z$avg_pi,pi.10kb.gu.z$avg_pi,pi.10kb.ru.w$avg_pi,pi.10kb.ty.w$avg_pi,pi.10kb.gu.w$avg_pi,outline=F)

### Scatterplot demo (for main figure)----------------------------------------

fstrho.100kb.ruty <- merge(fst.100kb.ruty.a,rho,by.x=c('chromosome','window_pos_1'),by.y=c('chrom','start.fix'))
dxyrho.100kb.ruty <- merge(dxy.100kb.ruty.a,rho,by.x=c('chromosome','window_pos_1'),by.y=c('chrom','start.fix'))
pirho.100kb.ru <- merge(pi.100kb.ru.a,rho,by.x=c('chromosome','window_pos_1'),by.y=c('chrom','start.fix'))
fstpi.100kb.ru <- merge(fst.100kb.ruty.a,pi.100kb.ru.a,by=c('chromosome','window_pos_1'))

par(mfrow=c(2,2))
plot(fstrho.100kb.ruty$rate,fstrho.100kb.ruty$avg_wc_fst,pch=20,col=alpha('#E28026',0.25),xlab='Recombination Rate',ylab='Fst')
plot(dxyrho.100kb.ruty$rate,dxyrho.100kb.ruty$avg_dxy,pch=20,col=alpha('#E28026',0.25),ylim=c(0,0.015),xlab='Recombination Rate',ylab='dxy')
plot(pirho.100kb.ru$rate,pirho.100kb.ru$avg_pi,pch=20,col=alpha('#B03160',0.25),ylim=c(0,0.015),xlab='Recombination Rate',ylab='pi')
plot(fstpi.100kb.ru$avg_pi,fstpi.100kb.ru$avg_wc_fst,pch=20,col=alpha('#B03160',0.25),xlim=c(0,0.015),xlab='pi',ylab='Fst')

### Scatterplots-------------------------------------------------------------

# Genome-wide-------
par(mfrow=c(3,3))
# π
plot(pi.100kb.ru$avg_pi,pi.100kb.ty$avg_pi,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(pi.100kb.ru$avg_pi,pi.100kb.gu$avg_pi,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(pi.100kb.ty$avg_pi,pi.100kb.gu$avg_pi,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
# dxy
plot(dxy.100kb.ruty$avg_dxy,dxy.100kb.rugu$avg_dxy,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(dxy.100kb.ruty$avg_dxy,dxy.100kb.guty$avg_dxy,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(dxy.100kb.rugu$avg_dxy,dxy.100kb.guty$avg_dxy,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
# Fst
plot(fst.100kb.ruty$avg_wc_fst,fst.100kb.rugu$avg_wc_fst,pch=20,col=alpha('black',0.15),xlim=c(0,0.9),ylim=c(0,0.9))
plot(fst.100kb.ruty$avg_wc_fst,fst.100kb.guty$avg_wc_fst,pch=20,col=alpha('black',0.15),xlim=c(0,0.9),ylim=c(0,0.9))
plot(fst.100kb.rugu$avg_wc_fst,fst.100kb.guty$avg_wc_fst,pch=20,col=alpha('black',0.15),xlim=c(0,0.9),ylim=c(0,0.9))

# π vs dxy
par(mfrow=c(2,3))
plot(pi.100kb.ru$avg_pi,dxy.100kb.ruty$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.ru$avg_pi,dxy.100kb.rugu$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.ty$avg_pi,dxy.100kb.ruty$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.ty$avg_pi,dxy.100kb.guty$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.gu$avg_pi,dxy.100kb.rugu$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.gu$avg_pi,dxy.100kb.guty$avg_dxy,pch=20,col=alpha('black',0.15))

# π vs Fst
pifst.100kb.ru.ruty <- merge(pi.100kb.ru,fst.100kb.ruty,by=c('chromosome','window_pos_1'))
pifst.100kb.ru.rugu <- merge(pi.100kb.ru,fst.100kb.rugu,by=c('chromosome','window_pos_1'))
pifst.100kb.ty.ruty <- merge(pi.100kb.ty,fst.100kb.ruty,by=c('chromosome','window_pos_1'))
pifst.100kb.ty.guty <- merge(pi.100kb.ty,fst.100kb.guty,by=c('chromosome','window_pos_1'))
pifst.100kb.gu.rugu <- merge(pi.100kb.gu,fst.100kb.rugu,by=c('chromosome','window_pos_1'))
pifst.100kb.gu.guty <- merge(pi.100kb.gu,fst.100kb.guty,by=c('chromosome','window_pos_1'))

plot(pifst.100kb.ru.ruty$avg_pi,pifst.100kb.ru.ruty$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.ru.rugu$avg_pi,pifst.100kb.ru.rugu$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.ty.ruty$avg_pi,pifst.100kb.ty.ruty$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.ty.guty$avg_pi,pifst.100kb.ty.guty$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.gu.rugu$avg_pi,pifst.100kb.gu.rugu$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.gu.guty$avg_pi,pifst.100kb.gu.guty$avg_wc_fst,pch=20,col=alpha('black',0.15))

# dxy vs Fst
dxyfst.100kb.ruty <- merge(dxy.100kb.ruty,fst.100kb.ruty,by=c('chromosome','window_pos_1'))
dxyfst.100kb.rugu <- merge(dxy.100kb.rugu,fst.100kb.rugu,by=c('chromosome','window_pos_1'))
dxyfst.100kb.guty <- merge(dxy.100kb.guty,fst.100kb.guty,by=c('chromosome','window_pos_1'))

plot(dxyfst.100kb.ruty$avg_dxy,dxyfst.100kb.ruty$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(dxyfst.100kb.rugu$avg_dxy,dxyfst.100kb.rugu$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(dxyfst.100kb.guty$avg_dxy,dxyfst.100kb.guty$avg_wc_fst,pch=20,col=alpha('black',0.15))

# rho vs Fst
rhofst.100kb.ruty <- merge(rho,fst.100kb.ruty,by.x=c('chrom','start.fix'),by.y=c('chromosome','window_pos_1'))
rhofst.100kb.rugu <- merge(rho,fst.100kb.rugu,by.x=c('chrom','start.fix'),by.y=c('chromosome','window_pos_1'))
rhofst.100kb.guty <- merge(rho,fst.100kb.guty,by.x=c('chrom','start.fix'),by.y=c('chromosome','window_pos_1'))

# Z chromosome-------
par(mfrow=c(3,3))
# π
plot(pi.100kb.ru.z$avg_pi,pi.100kb.ty.z$avg_pi,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(pi.100kb.ru.z$avg_pi,pi.100kb.gu.z$avg_pi,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(pi.100kb.ty.z$avg_pi,pi.100kb.gu.z$avg_pi,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
# dxy
plot(dxy.100kb.ruty.z$avg_dxy,dxy.100kb.rugu.z$avg_dxy,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(dxy.100kb.ruty.z$avg_dxy,dxy.100kb.guty.z$avg_dxy,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
plot(dxy.100kb.rugu.z$avg_dxy,dxy.100kb.guty.z$avg_dxy,pch=20,col=alpha('black',0.15),xlim=c(0,0.04),ylim=c(0,0.04))
# Fst
plot(fst.100kb.ruty.z$avg_wc_fst,fst.100kb.rugu.z$avg_wc_fst,pch=20,col=alpha('black',0.15),xlim=c(0,0.9),ylim=c(0,0.9))
plot(fst.100kb.ruty.z$avg_wc_fst,fst.100kb.guty.z$avg_wc_fst,pch=20,col=alpha('black',0.15),xlim=c(0,0.9),ylim=c(0,0.9))
plot(fst.100kb.rugu.z$avg_wc_fst,fst.100kb.guty.z$avg_wc_fst,pch=20,col=alpha('black',0.15),xlim=c(0,0.9),ylim=c(0,0.9))

# π vs dxy
par(mfrow=c(2,3))
plot(pi.100kb.ru.z$avg_pi,dxy.100kb.ruty.z$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.ru.z$avg_pi,dxy.100kb.rugu.z$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.ty.z$avg_pi,dxy.100kb.ruty.z$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.ty.z$avg_pi,dxy.100kb.guty.z$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.gu.z$avg_pi,dxy.100kb.rugu.z$avg_dxy,pch=20,col=alpha('black',0.15))
plot(pi.100kb.gu.z$avg_pi,dxy.100kb.guty.z$avg_dxy,pch=20,col=alpha('black',0.15))

# π vs Fst
pifst.100kb.ru.ruty.z <- merge(pi.100kb.ru.z,fst.100kb.ruty.z,by=c('chromosome','window_pos_1'))
pifst.100kb.ru.rugu.z <- merge(pi.100kb.ru.z,fst.100kb.rugu.z,by=c('chromosome','window_pos_1'))
pifst.100kb.ty.ruty.z <- merge(pi.100kb.ty.z,fst.100kb.ruty.z,by=c('chromosome','window_pos_1'))
pifst.100kb.ty.guty.z <- merge(pi.100kb.ty.z,fst.100kb.guty.z,by=c('chromosome','window_pos_1'))
pifst.100kb.gu.rugu.z <- merge(pi.100kb.gu.z,fst.100kb.rugu.z,by=c('chromosome','window_pos_1'))
pifst.100kb.gu.guty.z <- merge(pi.100kb.gu.z,fst.100kb.guty.z,by=c('chromosome','window_pos_1'))

plot(pifst.100kb.ru.ruty.z$avg_pi,pifst.100kb.ru.ruty.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.ru.rugu.z$avg_pi,pifst.100kb.ru.rugu.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.ty.ruty.z$avg_pi,pifst.100kb.ty.ruty.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.ty.guty.z$avg_pi,pifst.100kb.ty.guty.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.gu.rugu.z$avg_pi,pifst.100kb.gu.rugu.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(pifst.100kb.gu.guty.z$avg_pi,pifst.100kb.gu.guty.z$avg_wc_fst,pch=20,col=alpha('black',0.15))

# dxy vs Fst
dxyfst.100kb.ruty.z <- merge(dxy.100kb.ruty.z,fst.100kb.ruty.z,by=c('chromosome','window_pos_1'))
dxyfst.100kb.rugu.z <- merge(dxy.100kb.rugu.z,fst.100kb.rugu.z,by=c('chromosome','window_pos_1'))
dxyfst.100kb.guty.z <- merge(dxy.100kb.guty.z,fst.100kb.guty.z,by=c('chromosome','window_pos_1'))

plot(dxyfst.100kb.ruty.z$avg_dxy,dxyfst.100kb.ruty.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(dxyfst.100kb.rugu.z$avg_dxy,dxyfst.100kb.rugu.z$avg_wc_fst,pch=20,col=alpha('black',0.15))
plot(dxyfst.100kb.guty.z$avg_dxy,dxyfst.100kb.guty.z$avg_wc_fst,pch=20,col=alpha('black',0.15))

### Correlation coefficients-------------------------------------------------

# Genome-wide-------
# π
cor.test(pi.100kb.ru$avg_pi,pi.100kb.ty$avg_pi,method='spearman')
cor.test(pi.100kb.ru$avg_pi,pi.100kb.gu$avg_pi,method='spearman')
cor.test(pi.100kb.ty$avg_pi,pi.100kb.gu$avg_pi,method='spearman')
# dxy
cor.test(dxy.100kb.ruty$avg_dxy,dxy.100kb.rugu$avg_dxy,method='spearman')
cor.test(dxy.100kb.ruty$avg_dxy,dxy.100kb.guty$avg_dxy,method='spearman')
cor.test(dxy.100kb.rugu$avg_dxy,dxy.100kb.guty$avg_dxy,method='spearman')
# Fst
cor.test(fst.100kb.ruty$avg_wc_fst,fst.100kb.rugu$avg_wc_fst,method='spearman')
cor.test(fst.100kb.ruty$avg_wc_fst,fst.100kb.guty$avg_wc_fst,method='spearman')
cor.test(fst.100kb.rugu$avg_wc_fst,fst.100kb.guty$avg_wc_fst,method='spearman')

# π vs dxy
cor.test(pi.100kb.ru$avg_pi,dxy.100kb.ruty$avg_dxy,method='spearman')
cor.test(pi.100kb.ru$avg_pi,dxy.100kb.rugu$avg_dxy,method='spearman')
cor.test(pi.100kb.ty$avg_pi,dxy.100kb.ruty$avg_dxy,method='spearman')
cor.test(pi.100kb.ty$avg_pi,dxy.100kb.guty$avg_dxy,method='spearman')
cor.test(pi.100kb.gu$avg_pi,dxy.100kb.rugu$avg_dxy,method='spearman')
cor.test(pi.100kb.gu$avg_pi,dxy.100kb.guty$avg_dxy,method='spearman')

# π vs Fst
cor.test(pifst.100kb.ru.ruty$avg_pi,pifst.100kb.ru.ruty$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.ru.rugu$avg_pi,pifst.100kb.ru.rugu$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.ty.ruty$avg_pi,pifst.100kb.ty.ruty$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.ty.guty$avg_pi,pifst.100kb.ty.guty$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.gu.rugu$avg_pi,pifst.100kb.gu.rugu$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.gu.guty$avg_pi,pifst.100kb.gu.guty$avg_wc_fst,method='spearman')

# dxy vs Fst
cor.test(dxyfst.100kb.ruty$avg_dxy,dxyfst.100kb.ruty$avg_wc_fst,method='spearman')
cor.test(dxyfst.100kb.rugu$avg_dxy,dxyfst.100kb.rugu$avg_wc_fst,method='spearman')
cor.test(dxyfst.100kb.guty$avg_dxy,dxyfst.100kb.guty$avg_wc_fst,method='spearman')

# rho vs Fst
cor.test(as.numeric(rhofst.100kb.ruty$rate),rhofst.100kb.ruty$avg_wc_fst,method='spearman')
cor.test(as.numeric(rhofst.100kb.rugu$rate),rhofst.100kb.rugu$avg_wc_fst,method='spearman')
cor.test(as.numeric(rhofst.100kb.guty$rate),rhofst.100kb.guty$avg_wc_fst,method='spearman')


# Z chromosome-------
# π
cor.test(pi.100kb.ru.z$avg_pi,pi.100kb.ty.z$avg_pi,method='spearman')
cor.test(pi.100kb.ru.z$avg_pi,pi.100kb.gu.z$avg_pi,method='spearman')
cor.test(pi.100kb.ty.z$avg_pi,pi.100kb.gu.z$avg_pi,method='spearman')
# dxy
cor.test(dxy.100kb.ruty.z$avg_dxy,dxy.100kb.rugu.z$avg_dxy,method='spearman')
cor.test(dxy.100kb.ruty.z$avg_dxy,dxy.100kb.guty.z$avg_dxy,method='spearman')
cor.test(dxy.100kb.rugu.z$avg_dxy,dxy.100kb.guty.z$avg_dxy,method='spearman')
# Fst
cor.test(fst.100kb.ruty.z$avg_wc_fst,fst.100kb.rugu.z$avg_wc_fst,method='spearman')
cor.test(fst.100kb.ruty.z$avg_wc_fst,fst.100kb.guty.z$avg_wc_fst,method='spearman')
cor.test(fst.100kb.rugu.z$avg_wc_fst,fst.100kb.guty.z$avg_wc_fst,method='spearman')

# π vs dxy
cor.test(pi.100kb.ru.z$avg_pi,dxy.100kb.ruty.z$avg_dxy,method='spearman')
cor.test(pi.100kb.ru.z$avg_pi,dxy.100kb.rugu.z$avg_dxy,method='spearman')
cor.test(pi.100kb.ty.z$avg_pi,dxy.100kb.ruty.z$avg_dxy,method='spearman')
cor.test(pi.100kb.ty.z$avg_pi,dxy.100kb.guty.z$avg_dxy,method='spearman')
cor.test(pi.100kb.gu.z$avg_pi,dxy.100kb.rugu.z$avg_dxy,method='spearman')
cor.test(pi.100kb.gu.z$avg_pi,dxy.100kb.guty.z$avg_dxy,method='spearman')

# π vs Fst
cor.test(pifst.100kb.ru.ruty.z$avg_pi,pifst.100kb.ru.ruty.z$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.ru.rugu.z$avg_pi,pifst.100kb.ru.rugu.z$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.ty.ruty.z$avg_pi,pifst.100kb.ty.ruty.z$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.ty.guty.z$avg_pi,pifst.100kb.ty.guty.z$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.gu.rugu.z$avg_pi,pifst.100kb.gu.rugu.z$avg_wc_fst,method='spearman')
cor.test(pifst.100kb.gu.guty.z$avg_pi,pifst.100kb.gu.guty.z$avg_wc_fst,method='spearman')

# dxy vs Fst
cor.test(dxyfst.100kb.ruty.z$avg_dxy,dxyfst.100kb.ruty.z$avg_wc_fst,method='spearman')
cor.test(dxyfst.100kb.rugu.z$avg_dxy,dxyfst.100kb.rugu.z$avg_wc_fst,method='spearman')
cor.test(dxyfst.100kb.guty.z$avg_dxy,dxyfst.100kb.guty.z$avg_wc_fst,method='spearman')

### Summary statistics-------------------------------------------------------

# Genome-wide
mean(pi.100kb.ru$avg_pi,na.rm=T)
mean(pi.100kb.ty$avg_pi,na.rm=T)
mean(pi.100kb.gu$avg_pi,na.rm=T)
median(pi.100kb.ru$avg_pi,na.rm=T)
median(pi.100kb.ty$avg_pi,na.rm=T)
median(pi.100kb.gu$avg_pi,na.rm=T)
sd(pi.100kb.ru$avg_pi,na.rm=T)
sd(pi.100kb.ty$avg_pi,na.rm=T)
sd(pi.100kb.gu$avg_pi,na.rm=T)

mean(dxy.100kb.ruty$avg_dxy,na.rm=T)
mean(dxy.100kb.rugu$avg_dxy,na.rm=T)
mean(dxy.100kb.guty$avg_dxy,na.rm=T)
median(dxy.100kb.ruty$avg_dxy,na.rm=T)
median(dxy.100kb.rugu$avg_dxy,na.rm=T)
median(dxy.100kb.guty$avg_dxy,na.rm=T)
sd(dxy.100kb.ruty$avg_dxy,na.rm=T)
sd(dxy.100kb.rugu$avg_dxy,na.rm=T)
sd(dxy.100kb.guty$avg_dxy,na.rm=T)

mean(fst.100kb.ruty$avg_wc_fst,na.rm=T)
mean(fst.100kb.rugu$avg_wc_fst,na.rm=T)
mean(fst.100kb.guty$avg_wc_fst,na.rm=T)
median(fst.100kb.ruty$avg_wc_fst,na.rm=T)
median(fst.100kb.rugu$avg_wc_fst,na.rm=T)
median(fst.100kb.guty$avg_wc_fst,na.rm=T)
sd(fst.100kb.ruty$avg_wc_fst,na.rm=T)
sd(fst.100kb.rugu$avg_wc_fst,na.rm=T)
sd(fst.100kb.guty$avg_wc_fst,na.rm=T)

# Autosomes
mean(pi.100kb.ru.a$avg_pi,na.rm=T)
mean(pi.100kb.ty.a$avg_pi,na.rm=T)
mean(pi.100kb.gu.a$avg_pi,na.rm=T)
median(pi.100kb.ru.a$avg_pi,na.rm=T)
median(pi.100kb.ty.a$avg_pi,na.rm=T)
median(pi.100kb.gu.a$avg_pi,na.rm=T)
sd(pi.100kb.ru.a$avg_pi,na.rm=T)
sd(pi.100kb.ty.a$avg_pi,na.rm=T)
sd(pi.100kb.gu.a$avg_pi,na.rm=T)

mean(dxy.100kb.ruty.a$avg_dxy,na.rm=T)
mean(dxy.100kb.rugu.a$avg_dxy,na.rm=T)
mean(dxy.100kb.guty.a$avg_dxy,na.rm=T)
median(dxy.100kb.ruty.a$avg_dxy,na.rm=T)
median(dxy.100kb.rugu.a$avg_dxy,na.rm=T)
median(dxy.100kb.guty.a$avg_dxy,na.rm=T)
sd(dxy.100kb.ruty.a$avg_dxy,na.rm=T)
sd(dxy.100kb.rugu.a$avg_dxy,na.rm=T)
sd(dxy.100kb.guty.a$avg_dxy,na.rm=T)

mean(fst.100kb.ruty.a$avg_wc_fst,na.rm=T)
mean(fst.100kb.rugu.a$avg_wc_fst,na.rm=T)
mean(fst.100kb.guty.a$avg_wc_fst,na.rm=T)
median(fst.100kb.ruty.a$avg_wc_fst,na.rm=T)
median(fst.100kb.rugu.a$avg_wc_fst,na.rm=T)
median(fst.100kb.guty.a$avg_wc_fst,na.rm=T)
sd(fst.100kb.ruty.a$avg_wc_fst,na.rm=T)
sd(fst.100kb.rugu.a$avg_wc_fst,na.rm=T)
sd(fst.100kb.guty.a$avg_wc_fst,na.rm=T)

# Z chromosome
mean(pi.100kb.ru.z$avg_pi,na.rm=T)
mean(pi.100kb.ty.z$avg_pi,na.rm=T)
mean(pi.100kb.gu.z$avg_pi,na.rm=T)
median(pi.100kb.ru.z$avg_pi,na.rm=T)
median(pi.100kb.ty.z$avg_pi,na.rm=T)
median(pi.100kb.gu.z$avg_pi,na.rm=T)
sd(pi.100kb.ru.z$avg_pi,na.rm=T)
sd(pi.100kb.ty.z$avg_pi,na.rm=T)
sd(pi.100kb.gu.z$avg_pi,na.rm=T)

mean(dxy.100kb.ruty.z$avg_dxy,na.rm=T)
mean(dxy.100kb.rugu.z$avg_dxy,na.rm=T)
mean(dxy.100kb.guty.z$avg_dxy,na.rm=T)
median(dxy.100kb.ruty.z$avg_dxy,na.rm=T)
median(dxy.100kb.rugu.z$avg_dxy,na.rm=T)
median(dxy.100kb.guty.z$avg_dxy,na.rm=T)
sd(dxy.100kb.ruty.z$avg_dxy,na.rm=T)
sd(dxy.100kb.rugu.z$avg_dxy,na.rm=T)
sd(dxy.100kb.guty.z$avg_dxy,na.rm=T)

mean(fst.100kb.ruty.z$avg_wc_fst,na.rm=T)
mean(fst.100kb.rugu.z$avg_wc_fst,na.rm=T)
mean(fst.100kb.guty.z$avg_wc_fst,na.rm=T)
median(fst.100kb.ruty.z$avg_wc_fst,na.rm=T)
median(fst.100kb.rugu.z$avg_wc_fst,na.rm=T)
median(fst.100kb.guty.z$avg_wc_fst,na.rm=T)
sd(fst.100kb.ruty.z$avg_wc_fst,na.rm=T)
sd(fst.100kb.rugu.z$avg_wc_fst,na.rm=T)
sd(fst.100kb.guty.z$avg_wc_fst,na.rm=T)

# W chromosome
mean(pi.100kb.ru.w$avg_pi,na.rm=T)
mean(pi.100kb.ty.w$avg_pi,na.rm=T)
mean(pi.100kb.gu.w$avg_pi,na.rm=T)
median(pi.100kb.ru.w$avg_pi,na.rm=T)
median(pi.100kb.ty.w$avg_pi,na.rm=T)
median(pi.100kb.gu.w$avg_pi,na.rm=T)
sd(pi.100kb.ru.w$avg_pi,na.rm=T)
sd(pi.100kb.ty.w$avg_pi,na.rm=T)
sd(pi.100kb.gu.w$avg_pi,na.rm=T)

mean(dxy.100kb.ruty.w$avg_dxy,na.rm=T)
mean(dxy.100kb.rugu.w$avg_dxy,na.rm=T)
mean(dxy.100kb.guty.w$avg_dxy,na.rm=T)
median(dxy.100kb.ruty.w$avg_dxy,na.rm=T)
median(dxy.100kb.rugu.w$avg_dxy,na.rm=T)
median(dxy.100kb.guty.w$avg_dxy,na.rm=T)
sd(dxy.100kb.ruty.w$avg_dxy,na.rm=T)
sd(dxy.100kb.rugu.w$avg_dxy,na.rm=T)
sd(dxy.100kb.guty.w$avg_dxy,na.rm=T)

mean(fst.100kb.ruty.w$avg_wc_fst,na.rm=T)
mean(fst.100kb.rugu.w$avg_wc_fst,na.rm=T)
mean(fst.100kb.guty.w$avg_wc_fst,na.rm=T)
median(fst.100kb.ruty.w$avg_wc_fst,na.rm=T)
median(fst.100kb.rugu.w$avg_wc_fst,na.rm=T)
median(fst.100kb.guty.w$avg_wc_fst,na.rm=T)
sd(fst.100kb.ruty.w$avg_wc_fst,na.rm=T)
sd(fst.100kb.rugu.w$avg_wc_fst,na.rm=T)
sd(fst.100kb.guty.w$avg_wc_fst,na.rm=T)

