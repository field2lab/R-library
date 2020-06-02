#####################################convert to genlight##############################
library(LDcorSV)
library(adegenet)
library(ggplot2)
library("scales")
source("http://www.zzlab.net/GAPIT/emma.txt")
hapMap2genlight2 <- function(file){
  require(adegenet)
  hapmapraw1 <- read.csv(file.choose(), header=TRUE, row.names=1, sep="\t",
                         stringsAsFactors=FALSE,check.names = F)[,-(2:10)]
  monolist=c("A","C","G","T")
  delelist=c("-/+"," +/-")
  hapmapraw2 <- hapmapraw1[!hapmapraw1$alleles%in%monolist, ]
  hapmap <- rbind(hapmapraw2[!hapmapraw2$alleles%in%delelist, ],hapmapraw2[hapmapraw2$alleles%in%delelist, ])
  delecount=table(hapmapraw2$alleles%in%delelist)["TRUE"] 
  
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
  adele=s
  cdele=s
  gdele=s
  tdele=s
  dele=s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  names(adele) <- c("-","0","A","N")
  names(cdele) <- c("-","0","C","N")
  names(gdele) <- c("-","0","G","N")
  names(tdele) <- c("-","0","T","N")
  names(dele) <- c("-","0","Dummy","N")
  
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt,
               adele,adele,cdele,cdele,gdele,gdele,tdele,tdele,dele,dele)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G",
                   "A/-","-/A","C/-","-/C","G/-","-/G","T/-","-/T","-/+","+/-")
  
  # Pull out and convert genotypes
  S <- length(samples)
  SBlist <- vector(mode="list",S)   # set up list of SNPbin objects
  
  total <- S
  # create progress bar
  pb <- txtProgressBar(min = 1, max = S, style = 3)
  
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=T, USE.NAMES=FALSE)
    mygen[sapply(mygen , is.null)] <- NA
    # create SNPbin object for this individual
    SBlist[[i]] <- new("SNPbin", unlist(mygen))
    setTxtProgressBar(pb, i)
  }
  
  # make genlight object
  x <- new("genlight", SBlist,n.loc=length(loci),ind.names=samples,loc.names=loci)
  output=list(x,delecount)
  return(output)
}
mygeno=hapMap2genlight2(file.choose())
mygenodata=as.data.frame(mygeno[1])
hapmapinfo<- read.csv(file.choose(), header=TRUE, row.names=1, sep="\t",
                      stringsAsFactors=FALSE)[,1:3]
SNPlist=as.data.frame(as.matrix(unlist(colnames(mygenodata)),length(colnames(mygenodata)),3))
index <- match(SNPlist[,1], rownames(hapmapinfo))
SNPlist[,2]=hapmapinfo[index,2]
SNPlist[,3]=hapmapinfo[index,3]
colnames(SNPlist)=c("Marker","Chr","CM")
sum(!is.na(SNPlist[,2]))

LDdata=mygenodata
rownames(LDdata)=rownames(mygenodata)
LDmap=SNPlist[!is.na(SNPlist[,2]),]
colnames(LDmap)=c("Locus","LG","Position")
rownames(LDmap)=NULL
################LDdecay(limited distance)#################
LDall=NULL
for (i in 1:nlevels(as.factor(LDmap[,2]))){
LDdatatemp=LDdata[,colnames(LDdata)%in%LDmap[LDmap[,2]==levels(as.factor(LDmap[,2]))[i],1]]
LDmaptemp=LDmap[LDmap$LG==levels(as.factor(LDmap[,2]))[i],]
distancetemp=as.numeric(dist(LDmaptemp[,3]))
Locuslist=combn(LDmaptemp$Locus,2)
Locuslisttemp=Locuslist[,distancetemp<50000]

#foremmakinship=LDdata#use mygenodata if HWE is not an issue
#foremmakinship=rescale(as.matrix(foremmakinship), to=c(0,1))
#k_matrix=emma.kinship(t(foremmakinship),"additive","all") 
#rownames(k_matrix)=rownames(foremmakinship)
#colnames(k_matrix)=rownames(foremmakinship)

LDtemp=as.data.frame(NA,nrow=ncol(Locuslisttemp) ,ncol=3)
for (j in 1:ncol(Locuslisttemp)){
LDcal=cbind(LDdatatemp[,colnames(LDdatatemp)==Locuslisttemp[1,j]],LDdatatemp[,colnames(LDdatatemp)==Locuslisttemp[2,j]])
LDtemp[j,1]=Measure.R2(biloci = LDcal, na.presence=T)
LDtemp[j,2]=as.character(Locuslisttemp[1,j])
LDtemp[j,3]=as.character(Locuslisttemp[2,j])
}

index2 <- match(LDtemp[,2], SNPlist$Marker)
index3 <- match(LDtemp[,3], SNPlist$Marker)
LDtemp$distance=NA
LDtemp$distance=abs(as.numeric(SNPlist[index2,3])-as.numeric(SNPlist[index3,3]))
colnames(LDtemp)=c("r2","Locus1","Locus2","distance")
LDall=rbind(LDall,LDtemp)
LDall=LDall[order(LDall$r2,decreasing = T),]
rownames(LDtemp)=NULL
print(i)
}
################sampling distances#################
#
LDall500=LDall
#
LDall1000=LDall
#
LDall5000=LDall
#66
LDall10000=LDall
#175
LDall15000=LDall
#450
LDall20000=LDall
#1182
LDall30000=LDall
#1920
LDall40000=LDall

LDall50000=LDall
#4637
LDall70000=LDall
#
LDall200000=LDall


#66
lo110000 <-loess(LDall10000$r2~LDall10000$distance, span=0.9, degree=2)
pre10000 <- predict(lo110000)
frame10000 <- data.frame(LDall10000$distance, pre10000)
frame10000=frame10000[order(frame10000$LDall10000.distance),]
xval10000 <- frame10000[which.min(abs(0.2 - frame10000$pre10000)),]
#175
lo115000 <-loess(LDall15000$r2~LDall15000$distance, span=0.9, degree=2)
pre15000 <- predict(lo115000)
frame15000 <- data.frame(LDall15000$distance, pre15000)
frame15000=frame15000[order(frame15000$LDall15000.distance),]
xval15000 <- frame15000[which.min(abs(0.2 - frame15000$pre15000)),]
#450
lo120000 <-loess(LDall20000$r2~LDall20000$distance, span=0.9, degree=2)
pre20000 <- predict(lo120000)
frame20000 <- data.frame(LDall20000$distance, pre20000)
frame20000=frame20000[order(frame20000$LDall20000.distance),]
xval20000 <- frame20000[which.min(abs(0.2 - frame20000$pre20000)),]
#1182
lo130000 <-loess(LDall30000$r2~LDall30000$distance, span=0.9, degree=2)
pre30000 <- predict(lo130000)
frame30000 <- data.frame(LDall30000$distance, pre30000)
frame30000=frame30000[order(frame30000$LDall30000.distance),]
xval30000 <- frame30000[which.min(abs(0.2 - frame30000$pre30000)),]
#1920
lo140000 <-loess(LDall40000$r2~LDall40000$distance, span=0.9, degree=2)
pre40000 <- predict(lo140000)
frame40000 <- data.frame(LDall40000$distance, pre40000)
frame40000=frame40000[order(frame40000$LDall40000.distance),]
xval40000 <- frame40000[which.min(abs(0.2 - frame40000$pre40000)),]
#1920
lo150000 <-loess(LDall50000$r2~LDall50000$distance, span=0.9, degree=2)
pre50000 <- predict(lo150000)
frame50000 <- data.frame(LDall50000$distance, pre50000)
frame50000=frame50000[order(frame50000$LDall50000.distance),]
xval50000 <- frame50000[which.min(abs(0.2 - frame50000$pre50000)),]


png(file="LDtemp.png",width=12.57703333,height=8,units="in",res=300)
par(mar=c(5.1,5,4.1,2.1))
plot(1, type="n", xlab="distance (bp)", ylab="LD (r^2)", xlim=c(0, 50000), ylim=c(0, 0.6),cex.lab=2,
     cex.axis=2,
     cex.main=3)
points(LDall200000$distance,LDall200000$r2,pch = ".")
#lines(frame10000, col="blue")
#lines(frame15000, col="red")
#lines(frame20000, col="blue")
lines(frame30000, col="red",lwd=3)
lines(frame40000, col="blue",lwd=3)
lines(frame50000, col="green",lwd=3)
abline(h=0.2, col="black")

legend(30000, 0.55, legend=c("30000bp", "40000bp","50000bp"),title ="Maximum sample range",
       col=c("red", "blue","green"), lwd=2,cex=2,box.lty=0)


#abline(v=xval30000[1,1], col="orange")
#abline(v=xval40000[1,1], col="orange")
#abline(v=xval50000[1,1], col="orange")
#mtext(round(xval30000[1,1],2), side=1, line=0.05, at=xval30000[1,1], cex=0.9, col="orange")
#mtext(round(xval40000[1,1],2), side=1, line=0.05, at=xval40000[1,1], cex=0.9, col="orange")
#mtext(round(xval50000[1,1],2), side=1, line=0.05, at=xval50000[1,1], cex=0.9, col="orange")
dev.off()
round(xval30000[1,1],2)
round(xval40000[1,1],2)
round(xval50000[1,1],2)

p <- ggplot(LDtemp, aes(LDtemp$distance,LDtemp$r2)) + geom_point()

bin_size=5000
#bin_plot_dat <- function(bin_size){
nr_bins <- nrow(LDtemp) / bin_size
x2 <- rep(1:nr_bins * bin_size, each = bin_size)
y2 <- tapply(LDtemp$r2, x2, mean)
data.frame(x = unique(x2), y= y2)
#}

plot_dat2 <- bin_plot_dat(50)
p2 <- ggplot(plot_dat2, aes(x,y)) +
  geom_point()

p2 + geom_smooth()

#Fitting method2###########
# this data was obtained from genotyping of 242 individuals, so n = no. of gametes = no. of sampled chromosomes = 2 x 250
n = 2*242
## in my experience, changing the value of C does not do a whole lot of damage to your calculations
## you can try different values of Cs to see for yourself
## but reported C values range from about 0.5 to 2
Cstart <- c(C=1)
# let's fit a non linear model using the arbitrary C value
dist=LDall70000$distance
rsq=LDall70000$r2
modelC <- nls(r2 ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))), data=LDall70000, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter, 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ((10+rho*LDall70000$distance)/((2+rho*LDall70000$distance)*(11+rho*LDall70000$distance)))*(1+((3+rho*LDall70000$distance)*(12+12*rho*LDall70000$distance+(rho*LDall70000$distance)^2))/(n*(2+rho*LDall70000$distance)*(11+rho*LDall70000$distance)))

newfile1A <- data.frame(LDall70000$distance, newrsq)

#maxld <- max(file$rsq) #using max LD value from initial input file
maxld <- max(newfile1A$newrsq) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile1A$LDall70000.dist[which.min(abs(newfile1A$newrsq-halfdecay))]
newfile1A <- newfile1A[order(newfile1A$LDall70000.dist),]
f1 <- data.frame(newfile1A$LDall70000.dist, newfile1A$newrsq)
xval <- f1[which.min(abs(0.2 - f1$newfile1A.newrsq)),] #find x value where y=0.2
xval[1,1]
# create the plot
#png(file="LDall70000.png",width=12.57703333,height=8,units="in",res=700)
plot(LDall70000$distance, LDall70000$r2, pch=".", cex=2, xlab="distance (bp)", ylab="LD (r^2)")
lines(newfile1A$LDall70000.dist, newfile1A$newrsq, col="blue", lwd=2)
abline(h=0.2, col="red")
abline(v=xval[1,1], col="green")
mtext(round(xval[1,1],2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
#dev.off()

#######################LDblock(slide-window)######
LG=levels(as.factor(LDmap[,2]))
for (i in c(1,4,7,10,13,16,19)){
  LDdatatemp=LDdata[,colnames(LDdata)%in%LDmap[LDmap[,2]==levels(as.factor(LDmap[,2]))[i],1]]
  LDmaptemp=LDmap[LDmap$LG==levels(as.factor(LDmap[,2]))[i],]
  LD1LG=NULL
  for (j in 1:(round(nrow(LDmaptemp)/10)-1)){
distancetemp=LDmaptemp[j*10,]
LDdatarange=LDdatatemp[,((j-1)*10+1):((j+1)*10)]
LDcaltemp=LD.Measures(LDdatarange, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence=TRUE)
LDrangetemp=as.data.frame(cbind(distancetemp$LG,distancetemp$Position))
LDrangetemp$averageR2=mean(LDcaltemp$r2)
LD1LG=rbind(LD1LG,LDrangetemp)
print(paste(levels(as.factor(LDmap[,2]))[i],j,sep=" "))
  }
  colnames(LD1LG)=c("LG","position","average r2")
  LD1LG=LD1LG[order(LD1LG$`average r2`),]
  plot(LD1LG$position,LD1LG$`average r2`)
}