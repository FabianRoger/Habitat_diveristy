library(plyr)
library(ggplot2)

otutable<- read.delim("otutable.txt")
otutable$OTU_ID <- paste("OTU", otutable$OTU_ID,sep="_") 

ID <- read.delim("metadata.txt")
ID <- ID[,-6]

TAX <- read.table("OTU_TAX.txt", sep="\t", header=T, stringsAsFactors = F)
TAX <- TAX[,-9]
colnames(TAX)[1] <- "OTU_ID"
TAX$OTU_ID <- paste("OTU", TAX$OTU_ID,sep="_")

#merge TAX to otutable
OTUwT <- join(TAX,otutable)

# replace Bacteria in $Domain by Chloroplast for all Chloroplast sequences
OTUwT[OTUwT$Class=="Chloroplast",]$Domain<-"Chloroplast"

# replace Bacteria in $Domain by Mitochondria for all mitochondria sequences
OTUwT[OTUwT$Family=="mitochondria",]$Domain<-"Mitochondria"

# plot number of OTUs in each Domain
ggplot(OTUwT,aes(x=Domain,fill=Domain))+
  geom_bar(stat="bin")+
  labs(title="Distribution of OTUs")+
  theme_bw(base_size=15)


# plot number of Sequences in each Domain

#number of total reads
sum(OTUwT[,9:ncol(OTUwT)])

# read per domain
OTUprop<-ddply(OTUwT[,c(2,9:ncol(OTUwT))],.(Domain), function(x) sum(x[,-1]))

ggplot(OTUprop,aes(x=Domain,y=V1,fill=Domain))+
  geom_bar(stat="identity")+
  labs(y="total reads", title="Distribution of Reads")+
  theme_bw(base_size=15)

# remove all non-bacterial OTUs from OTUwT
OTUwToB<-OTUwT[OTUwT$Domain=="Bacteria",]

# number of discarded OTUs
nrow(otutable)-nrow(OTUwToB) # 322

# % of discarded OTUs
100-(nrow(OTUwToB)/nrow(otutable))*100 # 2.56 %

# which OTUs are discared? 
OTUd <- OTUwT[!OTUwT$OTU_ID %in% OTUwToB$OTU_ID, ]$OTU_ID

#export list of discarded OTUs
write.table(OTUd,"excluded_OTUs.txt",quote=F)

# number of discared sequences
sum(otutable[,-1])-sum(OTUwToB[,9:ncol(OTUwToB)]) #145944 ~ 6.5 %

write.table(OTUwToB,"OTU_onlyBAC_wTax.txt",sep="\t")

# remove all non-bacterial OTUs from OTU table without taxa information
otutable<-otutable[ otutable$OTU_ID %in% OTUwToB$OTU_ID, ]

#export clean OTU table without taxa information
write.table(otutable,"OTU_onlyBAC.txt",sep="\t")

#calculate read distribution among samples

Reads<-colSums(otutable[,-1])
Reads<-data.frame(sample_ID=names(Reads),ReadN=Reads)

# join metadate to Reads

Reads$sample_ID<-as.character(Reads$sample_ID)
Reads<-join(Reads,ID)

#plot depth by season and sediment
ggplot(Reads, aes(x=habitat, y=ReadN))+
  facet_wrap(~season)+
  geom_boxplot()

#plot depth with outliers removed
ggplot(Reads, aes(x=habitat, y=ReadN))+
  facet_wrap(~season)+
  scale_y_continuous(limits=c(0,60000))+
  geom_boxplot()

## no systematic bias for season or habitat

#plot by sample
ggplot(Reads, aes(x=hab_div,y=ReadN,colour=replicate))+
  facet_wrap(~habitat*season)+
  geom_point()+
  scale_y_continuous(limits=c(0,60000))+
  theme(legend.position="none")

