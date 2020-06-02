group.col=sample(colors()[-1], 20, replace = FALSE) 
group.col=c("yellowgreen","magenta4","blue","darkslategrey","lightseagreen","darkgoldenrod3","green3","cyan","burlywood4","darkred","orangered","purple","ivory4","plum4","springgreen2","hotpink4","black","lightskyblue4","darkorange","tan3","orchid4","sienna4")

hclustlist1=my3ddataset[my3ddataset[,16]==1,]
PCAclust1=PCAdata[my3ddataset[,16]==1,]
pca_matrixlist1 = scale(PCAclust1) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNAlist1 = apply(pca_matrixlist1, 2, function (x) {ifelse(is.na(x), 0, x)}) 
PCOADlist1=dist.gene(pca_matrix_noNAlist1 , method = "pairwise", pairwise.deletion = FALSE,
               variance = FALSE)
eigPCOAlist1=pcoa(PCOADlist1, correction="lingoes", rn=NULL)
PCOV1list1=eigPCOAlist1$vectors[,1]
PCOV2list1=eigPCOAlist1$vectors[,2]
PCOV3list1=eigPCOAlist1$vectors[,3]
biplot(pcoa(PCOADlist1, correction="none", rn=NULL), Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL)
myPCOAdatalist1=cbind(hclustlist1,PCOV1list1,PCOV2list1,PCOV3list1)
myPCOAdatalist1$Clade=as.factor(myPCOAdatalist1$Clade)
png(file="ECOdata.png",width=10,height=7,units="in",res=300)
ggplot(myPCOAdatalist1, aes(x=PCOV1list1, y=PCOV2list1,colour =State,group=State))+
  geom_point(size=3)+theme_bw()+xlab(paste("PC1(",pc1_varlist1,"%)",sep=""))+ylab(paste("PC2(",pc2_varlist1,"%)",sep=""))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"))+
  scale_colour_manual(name = "Clade",values = group.col)
dev.off()


hclustlist2=my3ddataset[my3ddataset[,16]==2,]
PCAclust2=PCAdata[my3ddataset[,16]==2,]
pca_matrixlist2 = scale(PCAclust2) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNAlist2 = apply(pca_matrixlist2, 2, function (x) {ifelse(is.na(x), 0, x)}) 
PCOADlist2=dist.gene(pca_matrix_noNAlist2 , method = "pairwise", pairwise.deletion = FALSE,
                     variance = FALSE)
eigPCOAlist2=pcoa(PCOADlist2, correction="lingoes", rn=NULL)
PCOV1list2=eigPCOAlist2$vectors[,1]
PCOV2list2=eigPCOAlist2$vectors[,2]
PCOV3list2=eigPCOAlist2$vectors[,3]
biplot(pcoa(PCOADlist2, correction="none", rn=NULL), Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL)
myPCOAdatalist2=cbind(hclustlist2,PCOV1list2,PCOV2list2,PCOV3list2)
myPCOAdatalist2$Clade=as.factor(myPCOAdatalist2$Clade)
png(file="ECOdata.png",width=10,height=7,units="in",res=300)
ggplot(myPCOAdatalist2, aes(x=PCOV1list2, y=PCOV2list2,colour =mapcode,group=mapcode))+
  geom_point(size=3)+theme_bw()+xlab(paste("PC1(",pc1_varlist2,"%)",sep=""))+ylab(paste("PC2(",pc2_varlist2,"%)",sep=""))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"))+
  scale_colour_manual(name = "Clade",values = group.col)
dev.off()


hclustlist3=my3ddataset[my3ddataset[,16]==3,]
PCAclust3=PCAdata[my3ddataset[,16]==3,]
pca_matrixlist3 = scale(PCAclust3) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNAlist3 = apply(pca_matrixlist3, 2, function (x) {ifelse(is.na(x), 0, x)}) 
PCOADlist3=dist.gene(pca_matrix_noNAlist3 , method = "pairwise", pairwise.deletion = FALSE,
                     variance = FALSE)
eigPCOAlist3=pcoa(PCOADlist3, correction="lingoes", rn=NULL)
PCOV1list3=eigPCOAlist3$vectors[,1]
PCOV2list3=eigPCOAlist3$vectors[,2]
PCOV3list3=eigPCOAlist3$vectors[,3]
biplot(pcoa(PCOADlist3, correction="none", rn=NULL), Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL)
myPCOAdatalist3=cbind(hclustlist3,PCOV1list3,PCOV2list3,PCOV3list3)
myPCOAdatalist3$Clade=as.factor(myPCOAdatalist3$Clade)
png(file="ECOdata.png",width=10,height=7,units="in",res=300)
ggplot(myPCOAdatalist3, aes(x=PCOV1list3, y=PCOV2list3,colour =State,group=State))+
  geom_point(size=3)+theme_bw()+xlab(paste("PC1(",pc1_varlist3,"%)",sep=""))+ylab(paste("PC2(",pc2_varlist3,"%)",sep=""))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16,face="bold"))+
  scale_colour_manual(name = "Clade",values = group.col)
dev.off()