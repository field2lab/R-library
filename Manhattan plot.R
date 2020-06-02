library(ggplot2)
library(qqman)

Data=read.csv(file.choose())
clab=c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D","UN")
mdata=Data[,1:4]
colnames(mdata)=c("SNP","CHR","BP","P")

png("Internode3_partitioning_2018.png", width=12, height=4, units="in", res=600)
par(bg=NA)
manhattan(mdata, col = c("blue4", "orange3"),suggestiveline = -log10(1e-03),genomewideline=F,chrlabs=clab,annotatePval=1e-03,annotateTop=T)
dev.off()

##dont worry about these#################
Data[Data$Chromosome==1,2]='1A' 
Data[Data$Chromosome==2,2]='1B' 
Data[Data$Chromosome==3,2]='1D' 
Data[Data$Chromosome==4,2]='2A' 
Data[Data$Chromosome==5,2]='2B' 
Data[Data$Chromosome==6,2]='2D'
Data[Data$Chromosome==7,2]='3A' 
Data[Data$Chromosome==8,2]='3B' 
Data[Data$Chromosome==9,2]='3D'
Data[Data$Chromosome==10,2]='4A' 
Data[Data$Chromosome==11,2]='4B' 
Data[Data$Chromosome==12,2]='4D'
Data[Data$Chromosome==13,2]='5A' 
Data[Data$Chromosome==14,2]='5B' 
Data[Data$Chromosome==15,2]='5D'
Data[Data$Chromosome==16,2]='6A' 
Data[Data$Chromosome==17,2]='6B' 
Data[Data$Chromosome==18,2]='6D'
Data[Data$Chromosome==19,2]='7A' 
Data[Data$Chromosome==20,2]='7B' 
Data[Data$Chromosome==21,2]='7D'