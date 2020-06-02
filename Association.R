library(Matrix)
library(lme4)
library(MASS)
library(lmerTest)

phenodata = read.csv(file.choose())
##########################################################################################################################################################
#boxcox
anova=aov(HD~ Populations,data=phenodata)
boxcox(anova,data=phenodata, lambda = seq(-10, 2, 1/10), 
       plotit = TRUE, eps = 1/50, xlab = expression(lambda),
       ylab = "log-Likelihood")
phenodata$newgplant=sqrt(phenodata$Plant.Wt.)
phenodata$newTN.no.plant=sqrt(phenodata$Tiller.No.)

##########################################################################################################################################################
#lsmeans
str(phenodata)
phenodata$Block=as.factor(phenodata$Block)
fit= lmer(DP1~ Populations+(1|Block),data=phenodata)
anova(fit)
lsmeans(fit)
lsmeans(fit, pairwise~"Variety.No.", adjust="tukey")
phenmeans=lsmeans(fit)$lsmeans.table
write.csv(phenmeans, file="phenlsmeans.csv")
new_phenodata = read.csv(file.choose())
pheno=new_phenodata[,3]
#########################################################################################################################################################
#imputation
table(is.na(genos_imputed))
genos=as.data.frame(PCAdata)
genos[is.na(genos)]="N"
genos_imputed=genos
for (i in 1:ncol(genos)){
  temp=genos[,i]
  zeroes=sum(temp==0,na.rm = T)
  ones=sum(temp==1,na.rm = T)
  twos=sum(temp==2,na.rm = T)
  avgones=ones/(zeroes+ones+twos)
  avgtwos=twos/(zeroes+ones+twos)
  avgzeroes=zeroes/(zeroes+ones+twos)
  if (avgones>avgtwos) {
    if(avgones>avgzeroes){genos_imputed[,i]=gsub("N",1,temp)}
    else{genos_imputed[,i]=gsub("N",0,temp)}
  }
  else{
    if(avgtwos>avgzeroes){genos_imputed[,i]=gsub("N",2,temp)}
    else{genos_imputed[,i]=gsub("N",0,temp)}
  }
  print(i)}
genos_imputed=apply(genos_imputed,2,as.numeric)
#########################################################################################################################################################
#raw pheno
ECOcleaned2A=read.table(file.choose(),header=T)
Famlist=read.table(file.choose())
rownames(Famlist)=Famlist[,2]
phenoraw= read.csv(file.choose())
mylist=read.csv(file.choose())
DAPClist=mylist[mylist[,2]%in%Famlist[,2],]

ECOcleaned2A=ECOcleaned2A[ECOcleaned2A[,2]%in%mylist[,2],]
newDAPClist=mylist[mylist[,2]%in%ECOcleaned2A[,2],]
associlist=newDAPClist
rownames(newDAPClist)=NULL
#associlist=DAPClist[!DAPClist[,2]=="PC22_101wrong",]

phenodata=phenoraw[phenoraw[,4]%in%associlist[,2],]
associSNP=genos_imputed
associSNP=associSNP[associlist[,2]%in%phenodata[,4],]
phenodata=phenodata[order(phenodata[,4]),]

phenotl=phenodata[!is.na(phenodata[,6]),][,6]
associSNPtl=associSNP[!is.na(phenodata[,6]),]

phenotw=phenodata[!is.na(phenodata[,7]),][,7]
associSNPtw=associSNP[!is.na(phenodata[,7]),]

phenohd=phenodata[!is.na(phenodata[,8]),][,8]
associSNPhd=associSNP[!is.na(phenodata[,8]),]

phenodp1=phenodata[!is.na(phenodata[,9]),][,9]
associSNPdp1=associSNP[!is.na(phenodata[,9]),]

phenodp2=phenodata[!is.na(phenodata[,10]),][,10]
associSNPdp2=associSNP[!is.na(phenodata[,10]),]

phenodp3=phenodata[!is.na(phenodata[,11]),][,11]
associSNPdp3=associSNP[!is.na(phenodata[,11]),]

phenoip1=phenodata[!is.na(phenodata[,12]),][,12]
associSNPip1=associSNP[!is.na(phenodata[,12]),]

phenoip2=phenodata[!is.na(phenodata[,13]),][,13]
associSNPip2=associSNP[!is.na(phenodata[,13]),]

phenoip3=phenodata[!is.na(phenodata[,14]),][,14]
associSNPip3=associSNP[!is.na(phenodata[,14]),]

phenoPL=associlist[!is.na(associlist[,10]),][,10]
associSNPPL=ECOcleaned2A[,7:ncol(ECOcleaned2A)]
#########################################################################################################################################################
pvals[which(pvals2==min(pvals2))]
#simple regression
pvals=c()
for (i in 1:ncol(associSNPPL)){
  pvals[i]=anova(lm(phenoPL~associSNPPL[,i]))[1,5]
  print(i)
}
plot(1:5033,-log10(pvals),type="p",cex.lab=0.8)
pvals[is.na(pvals)]=1
colnames(associSNPPL)[which(pvals==min(pvals))]
pvals[which(pvals==min(pvals))]
pvalsPL=matrix(unlist(pvals),nrow=5033,byrow=T)

colnames(associSNPPL[,which(pvals<0.000000000000000000000006)])
PloidylevelSNPs1=associSNPPL[,which(pvals<0.000000000000000000000006)]
rownames(PloidylevelSNPs1)=associlist[,2]
PloidylevelSNPs1$POP=associlist[,10]
PloidylevelSNPs1=rbind(PloidylevelSNPs1,pvals[which(pvals<0.000000000000000000000006)])
write.csv(PloidylevelSNPs1,file = "PloidySNP(simple regression).csv")


#simple regression + PCA
pca_matrix = scale(associSNPPL)
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})   
eig = prcomp(pca_matrix_noNA, center = T)
scores = eig$x
pc1=scores[,1]
pc2=scores[,2]
pc3=scores[,3]
pvals2=c()
for (i in 1:ncol(associSNPPL)){
  pvals2[i]=anova(lm(phenoPL~pc1+pc2+pc3+associSNPPL[,i]))[4,5]
  }
plot(1:5033,-log10(pvals2),type="p",cex.lab=0.8)
pvals2[is.na(pvals2)]=1
colnames(associSNPPL)[which(pvals2==min(pvals2))]
pvals2[which(pvals2==min(pvals2))]
colnames(associSNPPL[,which(pvals2<0.00040)])
PloidylevelSNPs2=associSNPPL[,which(pvals2<0.00040)]
rownames(PloidylevelSNPs2)=associlist[,2]
PloidylevelSNPs2$POP=associlist[,10]
PloidylevelSNPs2=rbind(PloidylevelSNPs2,pvals2[which(pvals2<0.00040)])
write.csv(PloidylevelSNPs2,file = "PloidySNP(simple regression+PCs).csv")

#logstic regression
phenoPL[phenoPL == 4]=0
phenoPL[phenoPL == 6]=1
phenoPL[phenoPL == 8]=1
as.data.frame(phenoPL)
pvals3=c()
for (i in 1:ncol(associSNPPL)){
  templogit=summary(glm(formula=phenoPL~associSNPPL[,i],family=binomial))
  ifelse(templogit$coefficients[2,4]=pvals3[i]=templogit$coefficients[2,4]
}
plot(1:402,-log10(pvals3),type="p")
dim(associSNPPL[,which(pvals3<0.00001)])
colnames(associSNPPL[,which(pvals3<0.0001)])
colnames(genos[,which(pvals3<0.0001)])
corploidy=cbind(associSNPPL[,which(pvals3<0.0001)],DAPClist[,10])
corploidyraw=cbind(genos[,which(pvals3<0.0001)],DAPClist[,10])
write.csv(corploidy,file = "Pacbioploidy.csv")

#logstic regression + PCA
pvals4=c()
for (i in 1:ncol(genos_imputed)){
  templogit=summary(glm(formula=pheno~pc1+pc2+pc3+genos_imputed[,i],family=binomial))
    pvals4[i]=templogit$coefficients[5,4]
 }
plot(1:402,-log10(pvals4),type="p")
dim(genos_imputed[,which(pvals4<0.01)])
colnames(genos_imputed[,which(pvals4<0.01)])
corploidy=cbind(genos_imputed[,which(pvals4<0.01)],DAPClist[,10])


#EMMA
source("/home/rorshach/Desktop/Workingg codes/emma.R") #loads all the functions in the emma.R file

k_matrix=emma.kinship(t(associSNPPL),"additive","all") #this will take a while
pvals3=emma.ML.LRT(t(phenohd),newgenos_imputed[1:536,],k_matrix)$ps 
plot(distance,-log10(pvals),,type="b", xlab="distance", ylab="Pvals")
points(distance,-log10(pvals2),type="b",col="green")
points(distance,-log10(pvals3),type="b",col="blue")
legend("bottomright",c("Models:","Naive","PCs","Emmer"),text.col=c("red","black","green","blue"),merge=TRUE)

#plot
plot(1:5033,-log10(pvals),type="p",ylim=c(0,15), cex.lab=1.5,cex.axis = 1.5)
points(1:5033,-log10(pvals2),type="p",col="blue",ylim=c(0,6), cex.lab=1.5,cex.axis = 1.5)


a1=-log10(runif(length(pvals),0,1))
a2=-log10(runif(length(pvals2),0,1))
a3=-log10(runif(length(pvals3),0,1))
b1=-log10(pvals)
b2=-log10(pvals2)
b3=-log10(pvals3)
qqplot(a1,b1,xlim=range(c(a1,a2)),ylim=range(c(b1,b2)),xlab="Expected Pvals",ylab="Models Pvals")
abline(0,1)
par(new=T)
qqplot(a2,b2,axes = FALSE,xlim=range(c(a1,a2)),ylim=range(c(b1,b2)),col="green",xlab=" ",ylab=" ")
par(new=T)
qqplot(a3,b3,axes = FALSE,xlim=range(c(a1,a3)),ylim=range(c(b1,b3)),col="blue",xlab=" ",ylab=" ")
legend("bottomright",c("Models:","Naive","PCA"),text.col=c("black","green"),merge=TRUE)

Contiglist_imputed=cbind(DAPClist[,2],DAPClist[,10],as.data.frame(genos_imputed[,which(pvals3<0.0001)]))
colnames(Contiglist_imputed)[1]='Population'
colnames(Contiglist_imputed)[2]='Ploidy'
Contiglist_NAs=cbind(DAPClist[,2],DAPClist[,10],genos[,which(pvals3<0.0001)])
colnames(Contiglist_NAs)[1]='Population'
colnames(Contiglist_NAs)[2]='Ploidy'
write.csv(Contiglist_imputed,file = "Contiglist_imputed.csv")
write.csv(Contiglist_NAs,file = "Contiglist_NAs.csv")