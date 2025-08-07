## July 30, 2014 - ADB
## This script obtains data for harp seal fecundity rates analysis
## It obtains it from multiple sources
## 1979 - 2007: Acces database
## 2008: Reconciliation between Access database and ovary xls files error checked by Lee
## 2009-2013:  ovary xls files error checked by Lee
## 2014: Reconciliation between ovary xls files error checked by Lee and field sampling data file

require(lubridate)
require(doBy)
require(reshape)                 

#### define cutoff data for early puppers
## 46: Feb 15  / 59: Feb 28
cutoff <- 51


###########################################################################
#################2014 Data################################
###########################################################################


### read field and lab (ov14) data and manipulate data frames
setwd('D:/Buren_files/MEGA/papersAle/harps_fecundity/data/')
field14 <- read.csv('2014/fieldsampling_2014.csv',header=T,as.is=T)
ov14 <- read.csv('2014/ovaries2014.csv',header=T,as.is=T)
field14 <- field14[,c(1:7,9:13)]  
hsf14 <- subset(field14, Species==1  & Sex=='F')

hsf14_preg <- subset(hsf14, Pregnant=="?")
hsf14$Mature[which(hsf14$Mature=="")] <-0
hsf14$Pregnant[which(hsf14$Pregnant=="")] <-0

names(hsf14)[names(hsf14)=="NewID"] <- "ID"

ov14 <- ov14[,c('ID','mat','emb')]
#ov14 <- ov14[,c('ID','mat','emb','aclx','alcax','iclx','ilcax')]

#merge field and lab data and remove everything prior to 2014
mat14 <- merge(hsf14,ov14,by='ID',all=T)
mat14 <- subset(mat14,ID >= 20140000)
mat14$Mature <- as.integer(mat14$Mature)

mat14$Pregnant[which(mat14$Pregnant=="?")] <- ifelse(mat14$emb[which(mat14$Pregnant=="?")]==1,1,0)
mat14$Pregnant <- as.integer(mat14$Pregnant)

## define maturity and pregnancy codes  ----  0=FALSE  --- 1=TRUE
## Note that columns mat and emb are from ovary xls files error checked by Lee 
## and columns Mature and Pregnant are from field sampling data file
mat14$maturity <-NA
mat14$maturity <- ifelse(mat14$Mature==0 | mat14$mat==1,0,NA)     
ind <- which(is.na(mat14$maturity))
mat14$maturity[ind] <- ifelse(mat14$Mature[ind]==1 | mat14$mat[ind]>1,1,10)
mat14$pregnancy <-NA
mat14$pregnancy <- ifelse(mat14$Pregnant==0 | mat14$emb==2,0,NA)     
ind <- which(is.na(mat14$pregnancy))
mat14$pregnancy[ind] <- ifelse(mat14$Pregnant[ind]==1 | mat14$emb[ind]==1,1,10)

## Add columns for consistency with next data frame and create data frame m14, which contains 2014 data. This will be used in conjunction  with data sets for other years
mat14$cohyear <- rep(2014,length(mat14$ID))
mat14$Age <- rep(NA,length(mat14$ID))
 m14 <- mat14[,c('ID','Species','Year','Month','Day','mat','emb','cohyear','maturity','pregnancy','Age')]

## remove all objects from workspace except data frame with 2014 data
rm(list=setdiff(ls(), c("m14","cutoff")))

###########################################################################
#################2008-2013 Data################################
###########################################################################
## Data from ovary and age xls files error checked by Lee
## name of age files
agefilenames <-   'XSEALA'
agefilenums <- c(1,3:12)


## read ovary data from csv file that Lee checked and combined for all years
ovary <- read.csv (file="D:/Buren_files/MEGA/papersAle/harps_fecundity/data/combined_ovary.csv", header=T)

#for (i in ovaryfilenums[2:length(ovaryfilenums)]) {
#     ov <- read.csv (file=paste(ovaryfilenames, formatC(i, width=2, flag="0"), ".csv",  sep=""), header=T)
#     ovary <- rbind(ovary,ov)}

## read age data from csv files that Lee checked
setwd('D:/Buren_files/MEGA/papersAle/harps_fecundity/data')     
age <-  read.csv (file=paste(agefilenames, formatC(1, width=2, flag="0"), ".csv",  sep=""), header=T)
for (i in agefilenums[2:length(agefilenums)]) {
     ag <- read.csv (file=paste(agefilenames, formatC(i, width=2, flag="0"), ".csv",  sep=""), header=T)
     age <- rbind(age,ag)}

## manipulate data sets     
age <- age[which(age$ID..>20079999),c('ID..','Age')]
names(age)[names(age)=="ID.."] <- "ID"                          
ovary <- ovary[which(ovary$Species==1),c('ID','Species','Year','Month','Day','mat','emb')]
ovary$cohyear <- as.integer(substr(ovary$ID,1,4))
ovary <- ovary[which(ovary$cohyear<2014 & ovary$cohyear>2007),]

## define maturity and pregnancy codes  ----  0=FALSE  --- 1=TRUE
ovary$maturity <-NA
ovary$maturity <- ifelse(ovary$mat==1,0,NA)
ind <- which(is.na(ovary$maturity))
ovary$maturity[ind] <- ifelse(ovary$mat[ind]>1,1,10)
ovary$pregnancy <-NA
ovary$pregnancy <- ifelse(ovary$emb==2,0,NA)
ind <- which(is.na(ovary$pregnancy))
ovary$pregnancy[ind] <- ifelse(ovary$emb[ind]==1,1,10)


## merge ovary and age data 2008-2013
ovaryage <- merge(ovary,age,by='ID',all.x=T)
## add 2014 data to 2008-2013
ovaryage <- rbind(ovaryage,m14)

## remove all objects from workspace except data frame with 2008-2014 data
rm(list=setdiff(ls(), c("ovaryage","cutoff")))


###########################################################################
#################1979-2008 Data################################
###########################################################################
## Data from Access database queries where I obtained ages for all harp seals and ovary data for all harp seals and exported to csv files

## read data
ages <- read.csv('D:/Buren_files/MEGA/papersAle/harps_fecundity/data/harps_ages.csv',header=T,as.is=T, strip.white=T)
fec <- read.csv('D:/Buren_files/MEGA/papersAle/harps_fecundity/data/harps_fecundity.csv',header=T,as.is=T, strip.white=T)

## manipulate data sets
ages <- ages[,setdiff(names(ages),'Species')]
fec$cohyear <- as.integer(substr(fec$ID,1,4))

## define maturity and pregnancy codes  ----  0=FALSE  --- 1=TRUE
fec$maturity <-NA
fec$maturity <- ifelse(fec$mat==1,0,NA)
ind <- which(is.na(fec$maturity))
fec$maturity[ind] <- ifelse(fec$mat[ind]>1,1,10)
fec$pregnancy <-NA
fec$pregnancy <- ifelse(fec$emb==2,0,NA)
ind <- which(is.na(fec$pregnancy))
fec$pregnancy[ind] <- ifelse(fec$emb[ind]==1,1,10)
ind <- which(fec$pregnancy==10)
fec$pregnancy[ind] <- ifelse(fec$mat[ind]==2 | fec$mat[ind]==3,1,0)


## merge age and ovary data sets
fec <- merge(fec,ages,all.x=T,by='ID')

## remove sex from ID
fec$ID <- as.integer(substr(fec$ID,1,8))

## find seals present in xls files (ovaryage object) and not present in database (fec object)
diff08 <- setdiff( ovaryage[(ovaryage$cohyear==2008),'ID'],  fec[(fec$cohyear==2008),'ID'])

## bind seals found in previous step to seals in database
fec <-rbind(fec, ovaryage[ovaryage$ID %in% diff08,])

## remove 2008 seals from ovaryage object (these are contained in fec object) and the  bind fec and ovaryage 
ovaryage <- ovaryage[which(ovaryage$cohyear>2008),]
ovaryage <- rbind(ovaryage,fec)

## remove all objects from workspace except data frame with 1979-2014 data
rm(list=setdiff(ls(), c("ovaryage","cutoff")))


###########################################################################
#################Data manipulation on entire data set################################
###########################################################################
## define early puppers
ovaryage$doy <- yday(with(ovaryage, as.Date(paste(formatC(Day, width=2, flag="0"),'/',formatC(Month, width=2, flag="0"),'/',Year,sep=''),format='%d/%m/%Y')))

ovaryage[which(ovaryage$doy>cutoff & ovaryage$doy<150 & ovaryage$mat==8),'pregnancy'] <- 1
ovaryage[which(ovaryage$doy>cutoff & ovaryage$doy<150 & ovaryage$mat==8),'maturity'] <- 1
ovaryage[which(ovaryage$doy>cutoff & ovaryage$doy<150 & ovaryage$mat==8),'mat'] <- 2

ovaryage$EP <- ifelse(ovaryage$mat==8,1,0)

## remove YOY and foetuses
# the 2nd line of this section removes all records based on Age: it removes all age 0s, '90' codes (i.e. foetus, starvling), and records that contain NAs in the column age
# the first line of this section stores the records that contain NAs in age in a separate object (as of today, there are multiple seals 2012-2014 that have not been aged)
# the last line puts the 2 objects together
oy <- ovaryage[which(is.na(ovaryage$Age)),]
ovaryage <- ovaryage[which(ovaryage$Age>0 & ovaryage$Age<90),]
ovaryage <- rbind(ovaryage,oy)

ovaryage <- subset(ovaryage, Month>9 | Month<3)

ovaryage <- ovaryage[order(ovaryage$ID),]
#write.csv(ovaryage, 'ovaryage.csv', row.names=F)

rm(list=setdiff(ls(), c("ovaryage","cutoff")))
###########################################################################
#################Data summaries################################
###########################################################################


matsum <-   summaryBy(maturity~cohyear+maturity,data=ovaryage,FUN=length)
pregsum <-   summaryBy(pregnancy~cohyear+pregnancy,data=ovaryage,FUN=length)
epsum <-   summaryBy(EP~cohyear+EP,data=ovaryage,FUN=length)

 meltmat <- melt(matsum, id=c("cohyear",'maturity','maturity.length'))
 mattable <- cast(meltmat,cohyear~maturity,value='maturity.length')
 mattable2 <- recast(matsum, cohyear~maturity,value='maturity.length', id.var=c("cohyear",'maturity','maturity.length'))
 names(mattable) <- c('cohyear','immature','mature')
 mattable$N <- with(mattable, mature+immature)
 
 
 meltpreg <- melt(pregsum, id=c("cohyear",'pregnancy','pregnancy.length'))
 pregtable <- cast(meltpreg,cohyear~pregnancy,value='pregnancy.length')
 names(pregtable) <- c('cohyear','nonpregnant','pregnant')
 pregtable <- pregtable[,c(1,3)]

 meltep <- melt(epsum, id=c("cohyear",'EP','EP.length'))
 eptable <- cast(meltep,cohyear~EP,value='EP.length')
 names(eptable) <- c('cohyear','nonEP','EP') 
 eptable <- eptable[,c(1,3)]

 fecdata <- merge(mattable,pregtable,by='cohyear')
 fecdata <- merge(fecdata,eptable,by='cohyear')
 fecdata <- fecdata[,c(1,4,3,2,5,6)]

#rm(list=setdiff(ls(), c("fecdata",'ovaryage')))

fecdata$fecrate<-fecdata$pregnant/fecdata$mature
fecdata$totpreg<-with(fecdata,EP+pregnant)
fecdata$abrate<-fecdata$EP/fecdata$totpreg

#print(fecdata)
cohyear <- data.frame(cohyear=1979:2014)
fecdata <- merge(fecdata,cohyear,all=T)

#write.csv(fecdata,'fecdata-apr22.csv',row.names=F,na='')


# write.csv(ovaryage,'fecdata_Stenson_etal_2015.csv',row.names=F,na='')

 #### Miscellaneous stuff - checking data
#
#nrow(ovar[which(ovar$cohyear==2009),])
#
#ov11 <- subset(ovaryage, cohyear>=2009 & cohyear<2011)
#ovar <-read.csv(infile,header=T,strip.white=T)
#ovar <- subset(ovar, Month>9 | Month<3)
#ovar <- ovar[,c('ID','maturity')]
#names(ovar) <- c('ID','matur')
#ovar$cohyear <- as.numeric(substr(ovar$ID, 1, 4))     
#ov11b <- subset(ovar, cohyear>=2009 & cohyear<2011)
#ov06b <- subset(ovar, cohyear==2006 & matur==8)
#ov06b$ID <- as.integer(substr(ov06b$ID, 1, 8))
#ov06 <- subset(ovaryage, cohyear==2006 & mat==8)
#compare <- merge(ov06, ov06b, by='ID',all=T)
#
#missale <- compare$ID[which(is.na(compare$cohyear.x))]
#
#
#dim(ov11b)
#compare <- merge(ov11, ov11b, by='ID',all=T)
#missdenis <- compare$ID[which(is.na(compare$cohyear.y))]
#missale <- compare$ID[which(is.na(compare$cohyear.x))]
#
##missale <- paste(missale, 'F',sep='')
#
#ma <- missale[1]
#for(i in 2:length(missale)){
#ma <- paste(ma,missale[i], sep=' Or ')}
#
#ovary <- read.csv (file="D:/Buren_files/MEGA/papersAle/harps_fecundity/data/combined_ovary.csv", header=T)
#ovary <- ovary[which(ovary$Species==1),c('ID','Species','Year','Month','Day','mat','emb')]
#ovary$cohyear <- as.integer(substr(ovary$ID,1,4))
#ovary <- ovary[which(ovary$cohyear<2014 & ovary$cohyear>2007),]
#
#
# ovar[which(ovar$ID %in% missale),]
# 
# 
# 
# ## name of morph files
#morphfilenames <-   'XSMD'
#morphfilenums <- 2:10
#
#setwd('D:/Buren_files/DFO/MM/Data/xlsx/Error Checked Files')     
#morph <-  read.csv (file=paste(morphfilenames, formatC(morphfilenums[1], width=2, flag="0"), ".csv",  sep=""), header=T)
#for (i in morphfilenums[2:length(morphfilenums)]) {
#     mo <- read.csv (file=paste(morphfilenames, formatC(i, width=2, flag="0"), ".csv",  sep=""), header=T)
#     morph <- rbind(morph,mo)}
#     
#morph <- morph[,c('ID','Pelage','Body.wt','Length')]
#
#
#morale <-  morph[which(morph$ID %in% missale),]
#ageale <-  age[which(age$ID %in% missale),]
#ale <- merge(morale,ageale,all=T)
#
#mordenis <-  morph[which(morph$ID %in% missdenis),]
#agedenis <-  age[which(age$ID %in% missdenis),]
#denis <- merge(mordenis,agedenis,all=T)
#
#
# ovaryage[which(ovaryage$ID %in% missdenis),]