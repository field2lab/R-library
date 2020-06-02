install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("EMMREML")
install.packages("scatterplot3d") 

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

GapitSNPdat=as.data.frame(read.csv(file.choose(), head = T))
Header=GapitSNPdat[,1:11]
Header[,3]=1
names(Header)=c("rs","alleles","chrom","pos","strand","assembly","center","protLSID","assayLSID","panel","QCcode")
colnames(GapitSNPdat)=gsub('-', '_', colnames(GapitSNPdat),fixed=T)
colnames(GapitSNPdat)=gsub('.', '_', colnames(GapitSNPdat),fixed=T)
colnames(GapitSNPdat)=gsub('X', '', colnames(GapitSNPdat),fixed=T)

GapitSNP=GapitSNPdat[,colnames(GapitSNPdat)%in%phenodata[,4]]

colnames(phenodata)
phenohd=phenodata[!is.na(phenodata[,8]),][,8]
Ghd=GapitSNP[,!is.na(phenodata[,8])]
Ghd1=cbind(Header,Ghd)
phenohd1=cbind(colnames(Ghd),as.data.frame(round(phenohd,digits = 3)))
colnames(phenohd1)=c("Taxa","hd")


myGAPIT <- GAPIT(
  Y=phenohd1,
  G=GapitSNPdat,
  PCA.total=3
)

str(phenohd1)
OutputSNP=cbind(Header,GapitSNP)

write.csv(OutputSNP,file = "OutputSNP.csv")
write.csv(phenodata,file = "Outputpheno.csv")


H=1-abs(PCAdata-1)
het.ind=as.data.frame(apply(H,1,mean,na.rm=T))
het.snp=as.data.frame(apply(H,2,mean,na.rm=T))
het.ind$clades=newDAPClist[,16]
het.ind$ploidy=newDAPClist[,10]
het.ind$names=newDAPClist[,9]
ylab.ind=paste("Frequency (out of ",length(het.ind)," individuals)",sep="")
ylab.snp=paste("Frequency (out of ",length(het.snp)," markers)",sep="")
pdf("GAPIT.Heterozygosity.pdf", width =10, height = 6)
par(mfrow=c(1,2),mar=c(5,5,1,1)+0.1)
box(het.ind,col="gray", main="",ylab=ylab.ind, xlab="Heterozygosity of individuals")
hist(het.snp,col="gray", main="",ylab=ylab.snp, xlab="Heterozygosity of markers")
dev.off()
par(mfrow=c(1,1))
par(cex.lab=2)

png(file="Heterozygosity_ploidy.png",width=8,height=7,units="in",res=500)
par(cex.axis=2.3)
boxplot(het.ind[!het.ind[,3]=='6',1] ~het.ind[!het.ind[,3]=='6',]$ploidy,boxlwd = 4,lwd=2, outline=F)
dev.off()

png(file="Heterozygosity_clade.png",width=8,height=7,units="in",res=500)
par(cex.axis=2.3)
boxplot(het.ind[,1] ~het.ind$clade,boxlwd = 4,lwd=2,outline=FALSE)
dev.off()
wilcox.test(het.ind[!het.ind[,3]=='6',1],het.ind[!het.ind[,3]=='6',]$ploidy,alternative="less" )
kruskal.test(het.ind[,1] ~ het.ind$clade, data = het.ind)
shapiro.test(het.ind[,1])

apply(H[185,],1, mean,na.rm = TRUE)
mean(H[185,])
str(H[185,])
