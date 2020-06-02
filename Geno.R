library(adegenet)
library(DMwR)
library(scrime)
library(ggplot2)
library(rgl)
library(dendextend)
hapMap2genlight <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)]
  samples <- names(hapmap)[-1]
  loci <- row.names(hapmap)
  
  # set up conversion table
  s <- as.integer(c(0,1,2,NA))
  ac <- s
  ag <- s
  at <- s
  cg <- s
  ct <- s
  gt <- s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G")
  
  # Pull out and convert genotypes
  S <- length(samples)
  SBlist <- vector(mode="list",S)   # set up list of SNPbin objects
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=TRUE, USE.NAMES=FALSE)
    # create SNPbin object for this individual
    SBlist[[i]] <- new("SNPbin", mygen)
  }
  
  # make genlight object
  x <- new("genlight", SBlist)
  locNames(x) <- loci
  indNames(x) <- samples
  
  return(x)
}
##########################################################################################################################################################
#input data
F1data <- hapMap2genlight(file.choose())
F1data2=as.matrix(mydata)

F1list=read.csv(file.choose())
rownames(F1data2)=F1list[,1]

for (i in 1:nrow(mydata2)){
  NaCountingforrow[i,1]=sum(is.na(mydata2[i,]))
  NaCountingforrow[i,2]=rownames(mydata2)[i]
}
write.csv(NaCountingforrow,file = "countingdata.csv")

str(NaCountingforrow)

table(mydata2[255,]==mydata2[254,])
#quality control
##########################################################################################################################################################
miss_ind=c()
for (i in 1:nrow(ECOdata)){
  miss_ind[i] = sum(is.na(ECOdata[i,]))/ncol(ECOdata)
}
hist(miss_ind)
newECOdata=ECOdata[miss_ind<0.2,]
miss_snp=c()
for (i in 1:ncol(newECOdata)){
  miss_snp[i] = sum(is.na(newECOdata[,i]))/ncol(newECOdata)
}
hist(miss_snp)
newECOdata=newECOdata[,miss_snp<0.05]

chisq_pval = c() 
for (i in 1:ncol(newECOdata)) {
  print(i)
  #Calculate the genotype frequencies for each of the three possible genotypes.  Note that "na.rm = T" must be used to exclude NAs
  homo_ref_sum = sum(newECOdata[,i] == 0, na.rm = T)
  hetero_sum = sum(newECOdata[,i] == 1, na.rm = T)
  homo_alt_sum = sum(newECOdata[,i] == 2,na.rm=T)
  #Get the sample size
  sample_size = sum(homo_ref_sum, hetero_sum, homo_alt_sum)
  #Calculate the allele frequencies for the two alleles.
  ref_freq = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
  alt_freq = ((homo_alt_sum * 2) + hetero_sum)/(sample_size * 2)
  #Calculate the expected genotype frequencies under Hardy-Weinberg Equilibrium.
  expected_homo_ref = ref_freq^2
  expected_hetero = 2*ref_freq*alt_freq
  expected_homo_alt = alt_freq^2
  #Perform a chisq.test to test for a difference between the expected and observed genotype frequencies.
  expected_geno_freqs = c(expected_homo_ref, expected_hetero, expected_homo_alt)
  observed_geno_freqs = c(homo_ref_sum, hetero_sum, homo_alt_sum)
  chisq_test = chisq.test(observed_geno_freqs, p = expected_geno_freqs)
  #the function chisq.test() returns a "list" of results. See help(chisq.test) under the heading "Value" for details of all results calculated by this function.
  #A "list" is a special type of object in R. In the command below, the "$" symbol is used to extract a p-value from a list of chisq.test results.
  chisq_pval[i] = chisq_test$p.value
}
hist (chisq_pval)
hist(chisq_pval,breaks=200)

afterHWEtest=newECOdata[,chisq_pval>=0.05]
newmylist=mylist[rownames(afterHWEtest),]

#Allele frequency test
frequency=c()
for (i in 1:ncol(mydata2)){
  frequency[i]=round(sum(mydata2[,i])/90,2)  
  print(i)
}
newmydata2=mydata2[,frequency>0.2&frequency<0.8]
#knnImputation(data, k = 5, scale = T, meth = "weighAvg", distData = NULL)
##########################################################################################################################################################
#PCA
pca_matrix = scale(mydata2) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})    #sets NAs to mean (0)
eig = prcomp(pca_matrix_noNA, center = T) #Performs PCA on pca_matrix_noNA
summary(eig) # print variance accounted for 
loadings(eig) # pc loadings 
plot(eig,type="lines") # scree plot 
screeplot(eig, type="lines")
eigenvalues = (eig$sdev)^2
scores = eig$x 
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
pc2_var = round((eigenvalues[2]/sum(eigenvalues)) * 100, 1)
pc3_var = round((eigenvalues[3]/sum(eigenvalues)) * 100, 1)
plot(scores[,1],scores[,2]) # Plot PC1 versus PC2
PCV1=scores[,1]
PCV2=scores[,2]
PCV3=scores[,3]
myPCAdata=cbind(mylistdata,PCV1,PCV2,PCV3)
rownames(myPCAdata)=NULL
#pedigree
################################################################################
myPCAdata$Crosses=factor(myPCAdata$Crosses, levels = c("IL 106 x IL 102", "PC 17-109 x IL 102", "IL 99 x IL 106", "IL 106 x IL 99", "PC 17-109 x IL 99","IL 106 x PC 17-109","PC 20-102 x PC 17-109"
                                                       ,"P1: IL102 107-1-2"
                                                       ,"P2: IL102 224 1-2"
                                                       ,"P3: IL102 224 1-2 B2"
                                                       ,"P4: IL102 322 1-2"
                                                       ,"P5: IL102 Original"
                                                       ,"P6: IL106 143-3-4"
                                                       ,"P7: IL99 139-2-4"
                                                       ,"P8: IL99 419-1-1"
                                                       ,"P9: PC17-109 145-1-2"
                                                       ,"P10: PC17-109 218-1-2"
                                                       ,"P11: PC20-102 223-2"))
png(file="PAG2016.png",width=10,height=7,units="in",res=300)
ggplot(myPCAdata, aes(x=PCV1, y=PCV2,shape=Crosses,colour =Crosses,group=Crosses))+
  geom_point(size=3)+theme_bw()+xlab("PC1(14.2%)")+ylab("PC2(7.5%)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"),legend.text = element_text(size = 10, hjust = 3, vjust = 3, face = 'bold'))+
  scale_colour_manual(name = "F1&Parents",values = c("darkslategray","green3","chocolate1","red","blue","yellow","purple",
                                                 "black","black","black","black","black","black","black","black","black","black","black"))+
  scale_shape_manual(name = "F1&Parents",values = c(16,16,16,16,16,16,16,2,3,4,5,6,7,8,9,0,11,12))
  # scale_shape_manual(name = "F1&Parents",values = c(18,18,18,18,18,18,18,25,22,4,8,22,0,24,"#","%",9,3))
dev.off() 
##########################################################################################################################################################
#imputation
newmydata2<- knncatimputeLarge(t(mydata2))
##########################################################################################################################################################
#second allele
myFreq <- glMean(mydata)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)
##########################################################################################################################################################
#first allele
myFreq <- glMean(mydata)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)
##########################################################################################################################################################
#missing allele
temp <- density(glNA(mydata), bw=10)
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
     xlim=c(0,1701))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(glNA(flu), rep(0, nLoc(flu)), pch="|", col="blue")
##########################################################################################################################################################
#missing allele for each family
missings=NA.posi(mydata)
##########################################################################################################################################################