library(rrBLUP)
library(glmnet)
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
mydata <- hapMap2genlight(file.choose())
mydata2=as.matrix(mydata)
head(mydata2)
#knnImputation&PCA(data, k = 5, scale = T, meth = "weighAvg", distData = NULL)
##########################################################################################################################################################
pca_matrix = scale(mydata2) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})    #sets NAs to mean (0)
markerset=pca_matrix_noNA
#Phenodataset
##########################################################################################################################################################
mylistdata=read.csv(file.choose(),head=T)
rownames(pca_matrix_noNA)=NULL
markerdata=pca_matrix_noNA
rownames(markerdata)=mylistdata[,c("Code")]
#genomic selection
##########################################################################################################################################################
BLUPsall2=BLUPsall[sample(nrow(BLUPsall)),]
#No.of folds
k=8
Fianltraitaccuracy=matrix(,nrow=6,ncol=2)
for (j in 1:6){
  traitaccuracy=matrix(,nrow=k)
for (i in 1:k){
  Wsize=round(nrow(BLUPsall2)/k-1)
  Pheno_valid=as.matrix(BLUPsall2[(1+Wsize*(i-1)):(Wsize*(i-1)+round(nrow(BLUPsall2)/k)),])
  train=setdiff(1:nrow(BLUPsall2),(1+Wsize*(i-1)):(Wsize*(i-1)+round(nrow(BLUPsall2)/k)))
  Pheno_train=BLUPsall2[train,]
  markerdata=markerdata[rownames(BLUPsall2),] 
  m_valid=markerdata[(1+Wsize*(i-1)):(Wsize*(i-1)+round(nrow(BLUPsall2)/k)),]
  m_train=markerdata[train,]

  trait=Pheno_train[,j]
  #rr-BLUP
  trait_M=mixed.solve(trait, Z=m_train, K=NULL,SE=F,return.Hinv=F)
  traitu=trait_M$u
  e=as.matrix(traitu)
  pred_trait_valid=m_valid%*%e
  pred_trait=(pred_trait_valid[,1])+trait_M$beta
  trait_valid=Pheno_valid[,j]
  trait_accuracy=cor(trait_valid,pred_trait_valid,use="complete")
  traitaccuracy[i,1]=trait_accuracy

}
Fianltraitaccuracy[j,1]=mean(traitaccuracy[,1])
Fianltraitaccuracy[j,2]=sd(traitaccuracy[,1]/sqrt(k))
}
rownames(Fianltraitaccuracy)=colnames(BLUPsall)
