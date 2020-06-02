##########################################################################################################################################################
#correlation matrix
library(Hmisc)
library(PerformanceAnalytics)
phenodata = read.csv("C:/Users/jiaguo3.UOFI/Box Sync/Box Sync/Jia/Joe's F1 nursery at Philo/Master.csv", header=T)
corrdata=phenodata[,13:21]
str(phenodata[,13:21])
cor(corrdata, use="na.or.complete",method='pearson')
mtrix=rcorr(as.matrix(corrdata), type="pearson")
write.table(mtrix$r,file="phenocorrelation matrix.csv",sep=",",row.names=F)

colnames(phenodata[,13:21]) =c("HD","Loding","Crown D.","Tiller L.","Stem D.","Plant Wt.","Tiller No.","Tiller Wt.","Leaf No.")
png(file="Phenomatrix1.png",width=10,height=7,units="in",res=500)
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
##########################################################################################################################################################
#boxplot
boxBLUPs=BLUPsall
boxBLUPs$Plant.Wt=ifelse(boxBLUPs[,1]<0,(-1)*sqrt(abs(boxBLUPs[,1])),sqrt(boxBLUPs[,1]))
newbox=boxBLUPs[,c(7,2,3,4,5,6)]
str(newbox)
png(file="Phenobox1.png",width=10,height=7,units="in",res=500)
par(font.axis = 2,cex.lab=2, cex.axis=1.3, mar=c(5,6,4,2)+0.1)
boxplot(newbox,col="darkcyan", ylab ="BLUP")
dev.off()
##########################################################################################################################################################
#PCA
variablenames <- c(names(as.data.frame(phenodata[,c(4,10,13:21)])))
str(phenodata)
means=aggregate(cbind(as.matrix(phenodata[,13:21]))~phenodata[,c("Code")]+phenodata[,c("family")], FUN = mean)
names(means) <-variablenames 
standardisedmeans <- as.data.frame(scale(means[3:11]))
rownames(standardisedmeans)=means[,1]

fit=princomp(standardisedmeans, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
screeplot(fit, type="lines")
prin1=fit$scores[,1]
prin2=fit$scores[,2]
prin3=fit$scores[,3]
prin4=fit$scores[,4]
prin5=fit$scores[,5]

PhenoPCAdata=cbind(means,prin1,prin2,prin3)

png("phenPCA.png", width=10, height=7, units="in", res=500)
ggplot(PhenoPCAdata, aes(x=prin1, y=prin2, colour =family, group=family))+
  geom_point(size=3)+theme_bw()+xlab("PC1(38.9%)")+ylab("PC2(23.2%)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"),legend.text = element_text( size =10, hjust = 3, vjust = 3, face = 'bold'))+
  scale_colour_manual(name = "F1 Families",values = c("darkslategray","green3","chocolate1","red","blue","yellow","purple"))
dev.off()

