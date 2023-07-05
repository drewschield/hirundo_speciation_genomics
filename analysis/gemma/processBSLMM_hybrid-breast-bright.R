#!/usr/bin/Rscript

## author: kira delmore, modified by evelien de greef

## Note: 03.21.23 - I'm trying to adopt this script from de Greef et al. 2023 on purple martins to combine results from multiple BSLMM chains

setwd('./output/')

# merging BSLMM outputs for analyses
id="gwas_full_bslmm.hybrid-breast-bright.run"

gemma.files <- list.files(pattern="gwas_full_bslmm.hybrid-breast-bright.run")
hyp.files<-grep(".hyp.txt",gemma.files,value=TRUE)
para.files<-grep(".param.txt",gemma.files,value=TRUE)

## look at run results
for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat.2 <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat.2 <- rbind(hyp.dat.2,dat)
  }
}

write.table(hyp.dat.2,paste0(id,"_hyp.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

pdf(paste0(id,"_hyp.pdf",sep=""))
hist(hyp.dat.2$h)
hist(hyp.dat.2$pve)
hist(hyp.dat.2$rho)
hist(hyp.dat.2$pge)
hist(hyp.dat.2$pi)
hist(hyp.dat.2$n_gam)
dev.off()

h     <- quantile(hyp.dat.2$h, probs = c(0.5, 0.025, 0.975)) # 50% quantile == median values, 0.025 and 0.975 are lower and upper 95%ETPI
PVE   <- quantile(hyp.dat.2$pve, probs = c(0.5, 0.025, 0.975))
rho   <- quantile(hyp.dat.2$rho, probs = c(0.5, 0.025, 0.975))
PGE   <- quantile(hyp.dat.2$pge, probs = c(0.5, 0.025, 0.975))
n_gam <- quantile(hyp.dat.2$n_gamma, probs = c(0.5,0.025, 0.975))

h_extra <- sd(hyp.dat.2$h)
PVE_extra <- sd(hyp.dat.2$pve)
rho_extra <- sd(hyp.dat.2$rho)
PGE_extra <- sd(hyp.dat.2$pge)
n_gam_extra <- sd(hyp.dat.2$n_gamma)
table_1 <- rbind(h_extra,PVE_extra,rho_extra,PGE_extra,n_gam_extra)

table_2 <- rbind(h,rho,PVE,PGE,n_gam)
table_3 <- cbind(table_2,table_1)
colnames(table_3) <- c("median", "low95EPTI", "upp95EPTI", "stdev")
write.table(table_3,paste0(id,"_summary_table.csv",sep=""),quote=FALSE,sep=",",row.names=TRUE,col.names=TRUE) ##could code this
rm(list=ls(pattern="dat"))

# Close big data
rm(hyp.dat.2)
rm(dat)

## look at SNP results
for (i in seq_along(para.files)) {
  if (i == 1) {
    para.dat <- read.table(para.files[i], header=TRUE)
  } else {
    dat <- read.table(para.files[i], header=TRUE)
    para.dat <- cbind(para.dat, dat[,4:7])
  }
}

para.mean <- data.frame(cbind(para.dat[,1:3],
                              apply(para.dat[,seq(5,dim(para.dat)[2],4)], MARGIN=1, FUN=mean),
                              apply(para.dat[,seq(6,dim(para.dat)[2],4)], MARGIN=1, FUN=mean),
                              apply(para.dat[,seq(7,dim(para.dat)[2],4)], MARGIN=1, FUN=mean)))

colnames(para.mean) <- c("CHR","RS","PS","alpha","beta","gamma")

para.mean$abs.beta.g  <- abs(para.mean$beta * para.mean$gamma) # calculate absolute beta x gamma
para.mean$dir.beta.g  <- para.mean$beta * para.mean$gamma # calculate directional beta x gamma

write.table(para.mean,paste0(id,"_para.txt",sep=""),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

para.mean$id=1:nrow(para.mean)

#pdf(paste0(id,"_para_BETA.pdf",sep=""))
#plot(para.mean$id,para.mean$abs.beta.g,col=as.numeric(para.mean$CHR),xlab="SNPs",ylab="BETA")
#dev.off()

#pdf(paste0(id,"_para_PIP.pdf",sep=""))
#plot(para.mean$id,para.mean$gamma,col=as.numeric(para.mean$CHR),xlab="SNPs",ylab="PIP")
#dev.off()

rm(list=ls(pattern="dat"))

sparse.effects <- para.mean[para.mean$abs.beta.g != 0,] # remove snps with no sparse effect on phenotypes
top1snp.effects   <- sparse.effects[which((sparse.effects$abs.beta.g > quantile(x=sparse.effects$abs.beta.g,probs=0.99))|(sparse.effects$dir.beta.g < quantile(x=sparse.effects$dir.beta.g,probs=0.01))),] # get top 1% high effect snps
top0.1snp.effects  <- sparse.effects[which((sparse.effects$abs.beta.g > quantile(x=sparse.effects$abs.beta.g,probs=0.999))|(sparse.effects$dir.beta.g < quantile(x=sparse.effects$dir.beta.g,probs=0.001))),] # get top 0.1% high effect snps
top0.01snp.effects  <- sparse.effects[which((sparse.effects$abs.beta.g > quantile(x=sparse.effects$abs.beta.g,probs=0.9999))|(sparse.effects$dir.beta.g < quantile(x=sparse.effects$dir.beta.g,probs=0.0001))),] #0.01% 

pip01       <- para.mean[which(para.mean$gamma >= 0.01),] # snps with gamma (i.e. PIP) > 0.01
pip10       <- para.mean[which(para.mean$gamma >= 0.10),] # gamma > 0.10
pip25       <- para.mean[which(para.mean$gamma >= 0.25),] # gamma > 0.25

write.table(pip01, paste0(id,"_PIP0.01.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
write.table(pip10, paste0(id,"_PIP0.10.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
write.table(pip25, paste0(id,"_PIP0.25.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
write.table(sparse.effects, paste0(id,"_sparse.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")

write.table(top1snp.effects, paste0(id,"_top1SNPs.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
write.table(top0.1snp.effects, paste0(id,"_top0.1SNPs.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
write.table(top0.01snp.effects, paste0(id,"_top0.01SNPs.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
