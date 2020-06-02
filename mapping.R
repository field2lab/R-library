# GAPIT - Genomic Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka, Feng Tian and You Tang
# Last update: September 15, 2015

#Install packages (Do this section only for new installation of R)
#-------------------------------------------------------------------------------
source("http://www.bioconductor.org/biocLite.R") 
biocLite("multtest")
#install.packages("gplots")
#install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d

#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#source("/Users/Zhiwu/Dropbox/Current/revolutionr/gapit/gapit_functions.txt")
#############################################################################################

is.even <- function(x) x %% 2 == 0
sumitY  <- read.csv("/Users/Rorshach/Box Sync/wheat/summit/GWAS/BLUP_C.csv", header=T, sep=",",
                    stringsAsFactors=FALSE,check.names = F)
sumitYlist  <- read.csv("/Users/Rorshach/Box Sync/wheat/summit/Previous/sumit/entry2genotype.csv", header=T, sep=",",
                        stringsAsFactors=FALSE,check.names = F)
sumitG <- read.csv("/Users/Rorshach/Box Sync/wheat/summit/Previous/genotype data/FLB2017AM_M0.05_TAXA0.5_SITE0.2.hmp.txt", header=F, sep="\t",
                   stringsAsFactors=FALSE)


colnames(sumitY)[1]="Taxa"
f1=function(x,y) x=y[x==y[,1],2]
sumitY[,1]=sapply(sumitY[,1], f1,y=sumitYlist,simplify = TRUE)
sumitY=sumitY[,-2]

sumitG[1,1:11]=c("rs","alleles","chrom", "pos", "strand", "assembly", "center", "protLSID", "assayLSID", "panel", "QCcode")
sumitG=sumitG[,sumitG[1,-(1:11)]%in%sumitY[,1]]

sumitY=sumitY[sumitY[,1]%in%sumitG[1,-(1:11)],]
target <- sumitG[1,-(1:11)]
sumitY=sumitY[match(target, sumitY$Taxa),]

#Run GAPIT
myGAPIT <- GAPIT(
  Y=sumitY,
  G=sumitG,
  PCA.total=3
)


myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)
``
#Using ECMLM by Li and et. al. (BMC Biology, 2014)
myGAPIT <- GAPIT(
  Y=sumitY,
  G=sumitG,
  PCA.total=3,
  kinship.cluster=c("average", "complete", "ward"),
  kinship.group=c("Mean", "Max"),
  group.from=200,
  group.to=1000000,
  group.by=10
)