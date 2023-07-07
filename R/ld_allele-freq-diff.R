############################################################################
# Barn swallow allele frequency differences between populations for LD scans
############################################################################

# To compare LD in candidate gene regions to a reasonable background, we need
# to look at sets of SNPs with matched allele frequency differences between
# parental populations (so that we are looking at similar slices of
# admixture LD). If there is coupling between candidate regions associated
# with mate choice traits, we would expect elevated LD between candidate
# regions compared to background SNPs with matched allele frequency differences.

# This script contains commands for parsing parental allele frequencies,
# calculating allele frequency differences, and parsing sets of candidate and
# background SNPs based on allele frequency differences.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./analysis/ld')

library(data.table)
library(tidyverse)

### Read in allele frequency data-------------------------------------------

af.ru <- read.table('./afd/allele-freq.all.rustica.txt',header=T)
af.ty <- read.table('./afd/allele-freq.all.tytleri.txt',header=T)
af.gu <- read.table('./afd/allele-freq.all.gutturalis.txt',header=T)

### Calculate allele frequency differences and set data tables--------------

afd.ruty <- c()
afd.ruty$CHROM <- af.ru$CHROM
afd.ruty$POS <- af.ru$POS
afd.ruty$AFD <- abs(af.ru$FREQ1-af.ty$FREQ1)
afd.ruty <- data.frame(afd.ruty)

afd.rugu <- c()
afd.rugu$CHROM <- af.ru$CHROM
afd.rugu$POS <- af.ru$POS
afd.rugu$AFD <- abs(af.ru$FREQ1-af.gu$FREQ1)
afd.rugu <- data.frame(afd.rugu)

afd.guty <- c()
afd.guty$CHROM <- af.ru$CHROM
afd.guty$POS <- af.ru$POS
afd.guty$AFD <- abs(af.ty$FREQ1-af.gu$FREQ1)
afd.guty <- data.frame(afd.guty)

### Write allele frequency difference tables---------------------------------

write.table(afd.ruty,'./afd/allele-freq-diff.all.rustica-tytleri.txt',quote=F,sep="\t",row.names=F)
write.table(afd.rugu,'./afd/allele-freq-diff.all.rustica-gutturalis.txt',quote=F,sep="\t",row.names=F)
write.table(afd.guty,'./afd/allele-freq-diff.all.tytleri-gutturalis.txt',quote=F,sep="\t",row.names=F)

# Close inputs:

rm(af.ru)
rm(af.ty)
rm(af.gu)
rm(afd.ruty)
rm(afd.rugu)
rm(afd.guty)

### Read in AFD in candidate and background regions--------------------------

# These inputs were parsed based on the AFD outputs above and the set of
# candidate trait associated regions using bedtools.

afd.ruty.cand <- read.table('./afd/allele-freq-diff.candidate.rustica-tytleri.txt',header=T)
afd.rugu.cand <- read.table('./afd/allele-freq-diff.candidate.rustica-gutturalis.txt',header=T)
afd.guty.cand <- read.table('./afd/allele-freq-diff.candidate.tytleri-gutturalis.txt',header=T)

# Try looking at significant SNPs in candidate regions only
#afd.ruty.cand <- read.table('./afd/allele-freq-diff.candidate.sig.rustica-tytleri.txt',header=T)

afd.ruty.back <- read.table('./afd/allele-freq-diff.background.rustica-tytleri.txt',header=T)
afd.rugu.back <- read.table('./afd/allele-freq-diff.background.rustica-gutturalis.txt',header=T)
afd.guty.back <- read.table('./afd/allele-freq-diff.background.tytleri-gutturalis.txt',header=T)

# Try looking at autosome-only backgrounds
#afd.ruty.back <- read.table('./afd/allele-freq-diff.background.autosome.rustica-tytleri.txt',header=T)
#afd.rugu.back <- read.table('./afd/allele-freq-diff.background.autosome.rustica-gutturalis.txt',header=T)
#afd.guty.back <- read.table('./afd/allele-freq-diff.background.autosome.tytleri-gutturalis.txt',header=T)

### Parse SNP sets by bins of AFD---------------------------------------------

# Candidate regions
afd.ruty.cand.2 <- afd.ruty.cand %>% filter(AFD >= 0.15 & AFD < 0.2)
afd.ruty.cand.25 <- afd.ruty.cand %>% filter(AFD >= 0.2 & AFD < 0.25)
afd.ruty.cand.3 <- afd.ruty.cand %>% filter(AFD >= 0.25 & AFD < 0.3)
afd.ruty.cand.35 <- afd.ruty.cand %>% filter(AFD >= 0.3 & AFD < 0.35)
afd.ruty.cand.4 <- afd.ruty.cand %>% filter(AFD >= 0.35 & AFD < 0.4)
afd.ruty.cand.45 <- afd.ruty.cand %>% filter(AFD >= 0.4 & AFD < 0.45)
afd.ruty.cand.5 <- afd.ruty.cand %>% filter(AFD >= 0.45 & AFD < 0.5)
afd.ruty.cand.55 <- afd.ruty.cand %>% filter(AFD >= 0.5 & AFD < 0.55)
afd.ruty.cand.6 <- afd.ruty.cand %>% filter(AFD >= 0.55 & AFD < 0.6)

afd.rugu.cand.2 <- afd.rugu.cand %>% filter(AFD >= 0.15 & AFD < 0.2)
afd.rugu.cand.25 <- afd.rugu.cand %>% filter(AFD >= 0.2 & AFD < 0.25)
afd.rugu.cand.3 <- afd.rugu.cand %>% filter(AFD >= 0.25 & AFD < 0.3)
afd.rugu.cand.35 <- afd.rugu.cand %>% filter(AFD >= 0.3 & AFD < 0.35)
afd.rugu.cand.4 <- afd.rugu.cand %>% filter(AFD >= 0.35 & AFD < 0.4)
afd.rugu.cand.45 <- afd.rugu.cand %>% filter(AFD >= 0.4 & AFD < 0.45)
afd.rugu.cand.5 <- afd.rugu.cand %>% filter(AFD >= 0.45 & AFD < 0.5)
afd.rugu.cand.55 <- afd.rugu.cand %>% filter(AFD >= 0.5 & AFD < 0.55)
afd.rugu.cand.6 <- afd.rugu.cand %>% filter(AFD >= 0.55 & AFD < 0.6)

afd.guty.cand.2 <- afd.guty.cand %>% filter(AFD >= 0.15 & AFD < 0.2)
afd.guty.cand.25 <- afd.guty.cand %>% filter(AFD >= 0.2 & AFD < 0.25)
afd.guty.cand.3 <- afd.guty.cand %>% filter(AFD >= 0.25 & AFD < 0.3)
afd.guty.cand.35 <- afd.guty.cand %>% filter(AFD >= 0.3 & AFD < 0.35)
afd.guty.cand.4 <- afd.guty.cand %>% filter(AFD >= 0.35 & AFD < 0.4)
afd.guty.cand.45 <- afd.guty.cand %>% filter(AFD >= 0.4 & AFD < 0.45)
afd.guty.cand.5 <- afd.guty.cand %>% filter(AFD >= 0.45 & AFD < 0.5)
afd.guty.cand.55 <- afd.guty.cand %>% filter(AFD >= 0.5 & AFD < 0.55)
afd.guty.cand.6 <- afd.guty.cand %>% filter(AFD >= 0.55 & AFD < 0.6)

# Backgrounds
afd.ruty.back.2 <- afd.ruty.back %>% filter(AFD >= 0.15 & AFD < 0.2)
afd.ruty.back.25 <- afd.ruty.back %>% filter(AFD >= 0.2 & AFD < 0.25)
afd.ruty.back.3 <- afd.ruty.back %>% filter(AFD >= 0.25 & AFD < 0.3)
afd.ruty.back.35 <- afd.ruty.back %>% filter(AFD >= 0.3 & AFD < 0.35)
afd.ruty.back.4 <- afd.ruty.back %>% filter(AFD >= 0.35 & AFD < 0.4)
afd.ruty.back.45 <- afd.ruty.back %>% filter(AFD >= 0.4 & AFD < 0.45)
afd.ruty.back.5 <- afd.ruty.back %>% filter(AFD >= 0.45 & AFD < 0.5)
afd.ruty.back.55 <- afd.ruty.back %>% filter(AFD >= 0.5 & AFD < 0.55)
afd.ruty.back.6 <- afd.ruty.back %>% filter(AFD >= 0.55 & AFD < 0.6)

afd.rugu.back.2 <- afd.rugu.back %>% filter(AFD >= 0.15 & AFD < 0.2)
afd.rugu.back.25 <- afd.rugu.back %>% filter(AFD >= 0.2 & AFD < 0.25)
afd.rugu.back.3 <- afd.rugu.back %>% filter(AFD >= 0.25 & AFD < 0.3)
afd.rugu.back.35 <- afd.rugu.back %>% filter(AFD >= 0.3 & AFD < 0.35)
afd.rugu.back.4 <- afd.rugu.back %>% filter(AFD >= 0.35 & AFD < 0.4)
afd.rugu.back.45 <- afd.rugu.back %>% filter(AFD >= 0.4 & AFD < 0.45)
afd.rugu.back.5 <- afd.rugu.back %>% filter(AFD >= 0.45 & AFD < 0.5)
afd.rugu.back.55 <- afd.rugu.back %>% filter(AFD >= 0.5 & AFD < 0.55)
afd.rugu.back.6 <- afd.rugu.back %>% filter(AFD >= 0.55 & AFD < 0.6)

afd.guty.back.2 <- afd.guty.back %>% filter(AFD >= 0.15 & AFD < 0.2)
afd.guty.back.25 <- afd.guty.back %>% filter(AFD >= 0.2 & AFD < 0.25)
afd.guty.back.3 <- afd.guty.back %>% filter(AFD >= 0.25 & AFD < 0.3)
afd.guty.back.35 <- afd.guty.back %>% filter(AFD >= 0.3 & AFD < 0.35)
afd.guty.back.4 <- afd.guty.back %>% filter(AFD >= 0.35 & AFD < 0.4)
afd.guty.back.45 <- afd.guty.back %>% filter(AFD >= 0.4 & AFD < 0.45)
afd.guty.back.5 <- afd.guty.back %>% filter(AFD >= 0.45 & AFD < 0.5)
afd.guty.back.55 <- afd.guty.back %>% filter(AFD >= 0.5 & AFD < 0.55)
afd.guty.back.6 <- afd.guty.back %>% filter(AFD >= 0.55 & AFD < 0.6)

# This results in way too many background SNPs for analysis;
# sample a random set equal to the number in candidate bins
afd.ruty.back.2 <- afd.ruty.back.2[sample(nrow(afd.ruty.back.2), length(afd.ruty.cand.2$AFD)),]
afd.ruty.back.25 <- afd.ruty.back.25[sample(nrow(afd.ruty.back.25), length(afd.ruty.cand.25$AFD)),]
afd.ruty.back.3 <- afd.ruty.back.3[sample(nrow(afd.ruty.back.3), length(afd.ruty.cand.3$AFD)),]
afd.ruty.back.35 <- afd.ruty.back.35[sample(nrow(afd.ruty.back.35), length(afd.ruty.cand.35$AFD)),]
afd.ruty.back.4 <- afd.ruty.back.4[sample(nrow(afd.ruty.back.4), length(afd.ruty.cand.4$AFD)),]
afd.ruty.back.45 <- afd.ruty.back.45[sample(nrow(afd.ruty.back.45), length(afd.ruty.cand.45$AFD)),]
afd.ruty.back.5 <- afd.ruty.back.5[sample(nrow(afd.ruty.back.5), length(afd.ruty.cand.5$AFD)),]
afd.ruty.back.55 <- afd.ruty.back.55[sample(nrow(afd.ruty.back.55), length(afd.ruty.cand.55$AFD)),]
afd.ruty.back.6 <- afd.ruty.back.6[sample(nrow(afd.ruty.back.6), length(afd.ruty.cand.6$AFD)),]

afd.rugu.back.2 <- afd.rugu.back.2[sample(nrow(afd.rugu.back.2), length(afd.rugu.cand.2$AFD)),]
afd.rugu.back.25 <- afd.rugu.back.25[sample(nrow(afd.rugu.back.25), length(afd.rugu.cand.25$AFD)),]
afd.rugu.back.3 <- afd.rugu.back.3[sample(nrow(afd.rugu.back.3), length(afd.rugu.cand.3$AFD)),]
afd.rugu.back.35 <- afd.rugu.back.35[sample(nrow(afd.rugu.back.35), length(afd.rugu.cand.35$AFD)),]
afd.rugu.back.4 <- afd.rugu.back.4[sample(nrow(afd.rugu.back.4), length(afd.rugu.cand.4$AFD)),]
afd.rugu.back.45 <- afd.rugu.back.45[sample(nrow(afd.rugu.back.45), length(afd.rugu.cand.45$AFD)),]
afd.rugu.back.5 <- afd.rugu.back.5[sample(nrow(afd.rugu.back.5), length(afd.rugu.cand.5$AFD)),]
afd.rugu.back.55 <- afd.rugu.back.55[sample(nrow(afd.rugu.back.55), length(afd.rugu.cand.55$AFD)),]
afd.rugu.back.6 <- afd.rugu.back.6[sample(nrow(afd.rugu.back.6), length(afd.rugu.cand.6$AFD)),]

afd.guty.back.2 <- afd.guty.back.2[sample(nrow(afd.guty.back.2), length(afd.guty.cand.2$AFD)),]
afd.guty.back.25 <- afd.guty.back.25[sample(nrow(afd.guty.back.25), length(afd.guty.cand.25$AFD)),]
afd.guty.back.3 <- afd.guty.back.3[sample(nrow(afd.guty.back.3), length(afd.guty.cand.3$AFD)),]
afd.guty.back.35 <- afd.guty.back.35[sample(nrow(afd.guty.back.35), length(afd.guty.cand.35$AFD)),]
afd.guty.back.4 <- afd.guty.back.4[sample(nrow(afd.guty.back.4), length(afd.guty.cand.4$AFD)),]
afd.guty.back.45 <- afd.guty.back.45[sample(nrow(afd.guty.back.45), length(afd.guty.cand.45$AFD)),]
afd.guty.back.5 <- afd.guty.back.5[sample(nrow(afd.guty.back.5), length(afd.guty.cand.5$AFD)),]
afd.guty.back.55 <- afd.guty.back.55[sample(nrow(afd.guty.back.55), length(afd.guty.cand.55$AFD)),]
afd.guty.back.6 <- afd.guty.back.6[sample(nrow(afd.guty.back.6), length(afd.guty.cand.6$AFD)),]

### Write candidate and background AFD matched outputs------------------------

write.table(afd.ruty.cand.2,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin2.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.25,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin25.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.3,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin3.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.35,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin35.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.4,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin4.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.45,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin45.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.5,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin5.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.55,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin55.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.cand.6,'./afd/allele-freq-diff.candidate.rustica-tytleri.bin6.txt',sep="\t",quote=F,row.names=F)

write.table(afd.rugu.cand.2,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin2.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.25,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin25.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.3,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin3.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.35,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin35.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.4,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin4.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.45,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin45.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.5,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin5.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.55,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin55.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.cand.6,'./afd/allele-freq-diff.candidate.rustica-gutturalis.bin6.txt',sep="\t",quote=F,row.names=F)

write.table(afd.guty.cand.2,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin2.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.25,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin25.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.3,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin3.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.35,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin35.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.4,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin4.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.45,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin45.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.5,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin5.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.55,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin55.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.cand.6,'./afd/allele-freq-diff.candidate.tytleri-gutturalis.bin6.txt',sep="\t",quote=F,row.names=F)

write.table(afd.ruty.back.2,'./afd/allele-freq-diff.background.rustica-tytleri.bin2.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.25,'./afd/allele-freq-diff.background.rustica-tytleri.bin25.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.3,'./afd/allele-freq-diff.background.rustica-tytleri.bin3.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.35,'./afd/allele-freq-diff.background.rustica-tytleri.bin35.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.4,'./afd/allele-freq-diff.background.rustica-tytleri.bin4.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.45,'./afd/allele-freq-diff.background.rustica-tytleri.bin45.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.5,'./afd/allele-freq-diff.background.rustica-tytleri.bin5.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.55,'./afd/allele-freq-diff.background.rustica-tytleri.bin55.txt',sep="\t",quote=F,row.names=F)
write.table(afd.ruty.back.6,'./afd/allele-freq-diff.background.rustica-tytleri.bin6.txt',sep="\t",quote=F,row.names=F)

write.table(afd.rugu.back.2,'./afd/allele-freq-diff.background.rustica-gutturalis.bin2.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.25,'./afd/allele-freq-diff.background.rustica-gutturalis.bin25.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.3,'./afd/allele-freq-diff.background.rustica-gutturalis.bin3.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.35,'./afd/allele-freq-diff.background.rustica-gutturalis.bin35.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.4,'./afd/allele-freq-diff.background.rustica-gutturalis.bin4.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.45,'./afd/allele-freq-diff.background.rustica-gutturalis.bin45.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.5,'./afd/allele-freq-diff.background.rustica-gutturalis.bin5.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.55,'./afd/allele-freq-diff.background.rustica-gutturalis.bin55.txt',sep="\t",quote=F,row.names=F)
write.table(afd.rugu.back.6,'./afd/allele-freq-diff.background.rustica-gutturalis.bin6.txt',sep="\t",quote=F,row.names=F)

write.table(afd.guty.back.2,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin2.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.25,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin25.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.3,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin3.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.35,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin35.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.4,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin4.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.45,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin45.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.5,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin5.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.55,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin55.txt',sep="\t",quote=F,row.names=F)
write.table(afd.guty.back.6,'./afd/allele-freq-diff.background.tytleri-gutturalis.bin6.txt',sep="\t",quote=F,row.names=F)
