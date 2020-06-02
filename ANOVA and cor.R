################ANOVA and lsmean###########################################################################################################################################
library(lmerTest)
library(pbkrtest)
library(car)
#options(contrasts=c("contr.sum","contr.poly"))
options(scipen=999)
#control=lmerControl(check.nlev.gtr.1 = "ignore")
is.even=function(x) x%%2==0
#import data
inputData <- read.csv(file.choose(),check.names = F);#Variable names won't change

# lets force Rto treat the variables as factors rather than integers
inputData$BLOCK<- factor(inputData$BLOCK)
inputData$ENTRY <- factor(inputData$ENTRY)
inputData$LOC <- factor(inputData$LOC)
inputData$YEAR <- factor(inputData$YEAR)
str(inputData)

# dataset combining two locations and two years
Estimates=matrix(NA,250,1)

#empty vector to store anova
Pvalues=c()
#loop for each variable
for (i in 7:19){
  #model (change the model if you have one year data or one location data)
  lm1 =lmer(inputData[,i]~ENTRY+(1|LOC)+(1|LOC:BLOCK)+(1|ENTRY:LOC)+(1|YEAR)+(1|ENTRY:YEAR)+(1|YEAR:LOC)+(1|YEAR:LOC:ENTRY),control = lmerControl(check.nlev.gtr.1 = "warning"), data=inputData)
  #model results
  k=lsmeansLT(lm1, test.effs="ENTRY")
  temp=cbind(as.character(k[[1]]$ENTRY),k[[1]]$Estimate,k[[1]][,3])
  #lsmeans
  Estimates=merge(Estimates, temp, by=1,all=TRUE)
  #Pvalues
  Pvalues=c(Pvalues,as.numeric(anova(lm1,ddfm="Kenward-Roger")[[6]]))
  print(i)
  #m.large  <- lmer(inputData[,i] ~ENTRY+(1|LOC)+(1|LOC:BLOCK)+(1|ENTRY:LOC)+(1|YEAR)+(1|ENTRY:YEAR)+(1|YEAR:LOC)+(1|YEAR:LOC:ENTRY), data=inputData)
  #m.small  <- lmer(inputData[,i] ~(1|LOC)+(1|LOC:BLOCK)+(1|ENTRY:LOC)+(1|YEAR)+(1|ENTRY:YEAR)+(1|YEAR:LOC)+(1|YEAR:LOC:ENTRY), data=inputData)
  #anova(m.large, m.small)
  #Anova(inputData.entry,type=c("III"), 
  #      test.statistic=c("F"))
  #KRmodcomp(m.large, m.small)
}
#yieldmean = as.data.frame(tapply(inputData$YIELD1_PL, inputData$ENTRY, na.rm=T, mean))
#plotdata=cbind(temp,yieldmean)
#plot(plotdata[,2], plotdata[,4] , col="blue")


#Calculate anova
colnames(Estimates)[c(is.even(1:ncol(Estimates)))]=c(colnames(inputData[7:22]))
write.csv(Estimates,file="combinedLSmean_SE.csv")
write.csv(Pvalues,file="combinedPvalues.csv")
################BLUP###################
BLUPEstimates=matrix(NA,250,1)
#loop for each variable
for (i in 7:19){
  #model (change the model if you have one year data or one location data)
  lm2 =lmer(inputData[,i] ~(1|ENTRY)+(1|LOC)+(1|LOC:BLOCK)+(1|ENTRY:LOC)+(1|YEAR)+(1|ENTRY:YEAR)+(1|YEAR:LOC)+(1|YEAR:LOC:ENTRY),control = lmerControl(check.nlev.gtr.1 = "warning"), data=inputData)
  #model results
  lmblup = ranef(lm2)
  temp=cbind(rownames(lmblup$ENTRY),lmblup$ENTRY)
  colnames(temp)[2]=colnames(inputData)[i]
  #lsmeans
  BLUPEstimates=merge(BLUPEstimates, temp, by=1,all=TRUE)
  print(i)
}
# save the brixlineblup output to a separate .csv file
write.csv(BLUPEstimates, file="SummitBLUPS.csv")

## Compare BLUP to line averages on a scatterplot
#temp=lmblup$ENTRY
#lmean = as.data.frame(tapply(inputData$ANTHESIS, inputData$ENTRY, na.rm=T, mean))
#plotdata=cbind(temp,lmean)
#plot(plotdata[,1], plotdata[,2] , col="blue")
################correlation matrix###########################################################################################################################################
library(Hmisc)
library(PerformanceAnalytics)
library("psych")
#import data
phenodata = read.csv(file.choose())
#subset data (Y1)
phenodataY1=phenodata[phenodata$YEAR==1,]
phenodataY1L1=phenodataY1[phenodataY1$LOC=="1",]
phenodataY1L2=phenodataY1[phenodataY1$LOC=="2",]
#subset data (Y2)
phenodataY2=phenodata[phenodata$YEAR==2,]
phenodataY2L1=phenodataY2[phenodataY2$LOC=="1",]
phenodataY2L2=phenodataY2[phenodataY2$LOC=="2",]
corrdata=phenodataY1L1[,c(11:43)]
#prepare dataset for rcorr
corrdata=phenodata[,c(11:20)]
#check variable property
str(corrdata)
#run correlation using corrdata
mtrix=rcorr(as.matrix(corrdata), type="pearson")
#write csv file to workspace
write.table(mtrix$r,file="summitcorrelation_selectcombined rmatrix.csv",sep=",",row.names=F)
write.table(mtrix$P,file="summitcorrelation_selectcombined pmatrix.csv",sep=",",row.names=F)
#nice graph (correlation matrix with distribution and r,p values)
png(file="Phenomatrix.png",width=10,height=7,units="in",res=500)
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="na.or.complete",method='pearson'))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  if(missing(cex.cor)) cex.cor <- 1.4/strwidth(txt)
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))
  color1=ifelse(round(test$p.value,3)<0.001,"red","black")
  color2=ifelse(r>0.4,"red","black")
  text(0.5, 0.25, paste("r=",txt),cex=1.7,col=color2)
  text(.5, .75, Signif,cex=1.7,col=color1)
}
panel.smooth<-function (x, y, col = "blue4", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 2))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="darkcyan", ...)
}

pairs(corrdata, upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = panel.smooth)
dev.off()
