setwd('D:/Buren_files/MEGA/papersAle/harps_fecundity/data/')
source('data_allyears.r')
library(gplots)
nboot <- 1000
ci <- 95
 memory.limit(size = 8000000)


oa <- ovaryage[,(colnames(ovaryage) %in% c('ID','cohyear','maturity','pregnancy','EP'))]
oa <- oa[order(oa$ID),]

## create matrices to store simulated values
bootab <- matrix(data=NA,ncol=nboot,nrow=length(unique(oa$cohyear)),dimnames=list(c(unique(oa$cohyear)),c(1:nboot)))
bootfec <- bootab
bootabpoint <- matrix(data=NA,nrow=length(unique(oa$cohyear)),dimnames=list(c(unique(oa$cohyear))))
bootfecpoint <- bootabpoint
i <- 1

for (y in unique(oa$cohyear)){
   for (nb in 1:nboot){
      dat <- subset(oa, cohyear==y)
      bootID <- data.frame(ID=sample(x=dat$ID, size=nrow(dat), replace = TRUE, prob = NULL))
      bootdat <- merge(dat, bootID, by='ID')
      bootmat <- aggregate(bootdat$maturity,by=list(bootdat$maturity),FUN=length)
      bootpreg <- aggregate(bootdat$pregnancy,by=list(bootdat$pregnancy),FUN=length)
      bootep <- aggregate(bootdat$EP,by=list(bootdat$EP),FUN=length)
      tmpfec <-  bootpreg[which(bootpreg$Group.1==1),'x']/bootmat[which(bootmat$Group.1==1),'x']
      bootfec[i,nb] <- ifelse(length(tmpfec)==0,0,tmpfec)
      tmpab <- bootep[which(bootep$Group.1==1),'x']/( bootep[which(bootep$Group.1==1),'x'] + bootpreg[which(bootpreg$Group.1==1),'x'])
      bootab[i,nb] <-  ifelse(length(tmpab)==0,0,tmpab)

      }


      bootmat <- aggregate(dat$maturity,by=list(dat$maturity),FUN=length)
      bootpreg <- aggregate(dat$pregnancy,by=list(dat$pregnancy),FUN=length)
      bootep <- aggregate(dat$EP,by=list(dat$EP),FUN=length)
      tmpfec <-  bootpreg[which(bootpreg$Group.1==1),'x']/bootmat[which(bootmat$Group.1==1),'x']
      bootfecpoint[i] <- ifelse(length(tmpfec)==0,0,tmpfec)
      tmpab <- bootep[which(bootep$Group.1==1),'x']/( bootep[which(bootep$Group.1==1),'x'] + bootpreg[which(bootpreg$Group.1==1),'x'])
      bootabpoint[i] <-  ifelse(length(tmpab)==0,0,tmpab)

    i <- i +1
   }

feclim <- as.data.frame(t(apply(bootfec, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE)))
names(feclim) <- c('feclb','fecub')
feclim$cohyear <-   unique(oa$cohyear)

ablim <- as.data.frame(t(apply(bootab, 1, quantile, probs = c((1-ci/100)/2, ci/100+(1-ci/100)/2),  na.rm = TRUE)))
names(ablim) <- c('ablb','abub')
ablim$cohyear <-   unique(oa$cohyear)



bootCI <- data.frame(feclim,ablim)#
bootCI <- bootCI[,c(3,1,2,4,5)]
#write.csv(bootCI ,'bootCI.csv',row.names=F)
#par(mfrow=c(1,2))
#with(bootCI,plotCI(x=cohyear,y=abmean,ui=abub,li=ablb,type='o',lty=1,sfrac=0,pch=16,gap=0))
#with(bootCI,plotCI(x=cohyear,y=fecmean,ui=fecub,li=feclb,type='o',lty=1,sfrac=0,pch=16,gap=0))
#
#pointest <- data.frame(cohort_year=unique(oa$cohyear),fecrate=bootfecpoint,abrate= bootabpoint)
