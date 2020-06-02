library(randomcoloR)
brewer.pal(7, "BrBG")
group.col=sample(colors()[-1], 14, replace = FALSE) 
group.col=c("yellowgreen","magenta4","blue","darkslategrey","lightseagreen","darkgoldenrod3","green3","cyan","burlywood4","darkred","orangered","purple","ivory4","plum4","springgreen2","hotpink4","black","lightskyblue4","darkorange","tan3","orchid4","sienna4")

my3ddataset[,8]=as.character(my3ddataset[,8])
my3ddataset[,10]=as.character(my3ddataset[,10])

my3ddataset[88,8]='KS'
my3ddataset[89,8]='KS'
my3ddataset[90,8]='KS'
my3ddataset[91,8]='KS'
my3ddataset[92,8]='KS'

my3ddataset[108,8]='Red River'
my3ddataset[109,8]='Red River'


hclustlist1=my3ddataset[my3ddataset[,16]==1,]
PCAclust1=PCAdata[my3ddataset[,16]==1,]
pca_matrixlist1 = scale(PCAclust1) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNAlist1 = apply(pca_matrixlist1, 2, function (x) {ifelse(is.na(x), 0, x)}) 
eiglist1 = prcomp(pca_matrix_noNAlist1, center = T) #Performs PCA on pca_matrix_noNA
summary(eiglist1) # print variance accounted for 
loadings(eiglist1) # pc loadings 
plot(eiglist1,type="lines") # scree plot 
eigenvalueslist1 = (eiglist1$sdev)^2
pc1_varlist1 = round((eigenvalueslist1[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
pc2_varlist1 = round((eigenvalueslist1[2]/sum(eigenvalues)) * 100, 1)
pc3_varlist1= round((eigenvalueslist1[3]/sum(eigenvalues)) * 100, 1)
plot(eiglist1$x[,1],eiglist1$x[,2]) # Plot PC1 versus PC2
PCV1list1=eiglist1$x[,1]
PCV2list1=eiglist1$x[,2]
PCV3list1=eiglist1$x[,3]
myPCAdatalist1=cbind(hclustlist1,PCV1list1,PCV2list1,PCV3list1)
rownames(myPCAdatalist1)=NULL
myPCAdatalist1$Clade=as.factor(myPCAdatalist1$Clade)
g1_1 <- subset(myPCAdatalist1, State2== "NY1")
g1_2 <- subset(myPCAdatalist1, State2== "NY2")
myPCAdatalist1=myPCAdatalist1[-c(76:79),]
png(file="Clade1.png",width=10,height=7,units="in",res=300)
ggplot(myPCAdatalist1, aes(x=PCV1list1, y=PCV2list1,shape=Ploidy,colour =State,group=State))+
  geom_point(size=6)+theme_bw()+xlab(paste("PC1(",pc1_varlist1,"%)",sep=""))+ylab(paste("PC2(",pc2_varlist1,"%)",sep=""))+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title=element_text(size=16,face="bold"))+
  geom_point(data=g1_1[1,], colour="Black", shape=c(7), size=6) +
  geom_point(data=g1_2, colour="Black", shape=c(9), size=6) +
  theme(legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = 0)+
  theme(legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))+
  scale_x_continuous(limits = c(-45, 70))+
  scale_y_continuous(limits = c(-65, 50))+
  scale_colour_manual(name = "Collection site",values = group.col)+
scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(color = guide_legend(order=2),
          shape = guide_legend(order=1))
dev.off()
'c(6) is STP
c(2) is KST'
'scale_fill_manual(values = group.col,guide = guide_legend(reverse = TRUE))'

hclustlist2=my3ddataset[my3ddataset[,16]==2,]
PCAclust2=PCAdata[my3ddataset[,16]==2,]
pca_matrixlist2 = scale(PCAclust2) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNAlist2 = apply(pca_matrixlist2, 2, function (x) {ifelse(is.na(x), 0, x)}) 
eiglist2 = prcomp(pca_matrix_noNAlist2, center = T) #Performs PCA on pca_matrix_noNA
summary(eiglist2) # print variance accounted for 
loadings(eiglist2) # pc loadings 
plot(eiglist2,type="lines") # scree plot 
eigenvalueslist2 = (eiglist2$sdev)^2
pc1_varlist2 = round((eigenvalueslist2[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
pc2_varlist2 = round((eigenvalueslist2[2]/sum(eigenvalues)) * 100, 1)
pc3_varlist2= round((eigenvalueslist2[3]/sum(eigenvalues)) * 100, 1)
plot(eiglist2$x[,1],eiglist2$x[,2]) # Plot PC1 versus PC2
PCV1list2=eiglist2$x[,1]
PCV2list2=eiglist2$x[,2]
PCV3list2=eiglist2$x[,3]
myPCAdatalist2=cbind(hclustlist2,PCV1list2,PCV2list2,PCV3list2)
rownames(myPCAdatalist2)=NULL
myPCAdatalist2$Clade=as.factor(myPCAdatalist2$Clade)
png(file="Clade2.png",width=10,height=7,units="in",res=300)
ggplot(myPCAdatalist2, aes(x=PCV1list2, y=PCV2list2,colour =State,group=State,shape=Ploidy))+
  geom_point(size=6)+theme_bw()+xlab(paste("PC1(",pc1_varlist2,"%)",sep=""))+ylab(paste("PC2(",pc2_varlist2,"%)",sep=""))+
  theme(axis.text=element_text(size=10),        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = 0)+
  theme(legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))+
  scale_x_continuous(limits = c(-45, 70))+
  scale_y_continuous(limits = c(-65, 50))+
  scale_colour_manual(name = "Collection site",values = group.col)+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1))
dev.off()

hclustlist3=my3ddataset[my3ddataset[,16]==3,]
PCAclust3=PCAdata[my3ddataset[,16]==3,]
pca_matrixlist3 = scale(PCAclust3) #adjusts mean and stdev of each column to 1 and 0 respectively
pca_matrix_noNAlist3 = apply(pca_matrixlist3, 2, function (x) {ifelse(is.na(x), 0, x)}) 
eiglist3 = prcomp(pca_matrix_noNAlist3, center = T) #Performs PCA on pca_matrix_noNA
summary(eiglist3) # print variance accounted for 
loadings(eiglist3) # pc loadings 
plot(eiglist3,type="lines") # scree plot 
eigenvalueslist3 = (eiglist3$sdev)^2
pc1_varlist3 = round((eigenvalueslist3[1]/sum(eigenvalues)) * 100, 1) #% variance explained by PC1
pc2_varlist3 = round((eigenvalueslist3[2]/sum(eigenvalues)) * 100, 1)
pc3_varlist3= round((eigenvalueslist3[3]/sum(eigenvalues)) * 100, 1)
plot(eiglist3$x[,1],eiglist3$x[,2]) # Plot PC1 versus PC2
PCV1list3=eiglist3$x[,1]
PCV2list3=eiglist3$x[,2]
PCV3list3=eiglist3$x[,3]
myPCAdatalist3=cbind(hclustlist3,PCV1list3,PCV2list3,PCV3list3)
rownames(myPCAdatalist3)=NULL
myPCAdatalist3$Clade=as.factor(myPCAdatalist3$Clade)
g3 <- subset(myPCAdatalist3, State == "Red River")
myPCAdatalist3=myPCAdatalist3[-c(31,32),]
png(file="Clade3.png",width=10,height=7,units="in",res=300)
ggplot(myPCAdatalist3, aes(x=PCV1list3, y=PCV2list3,colour =State,group=State, shape=Ploidy))+
  geom_point(size=6)+theme_bw()+xlab(paste("PC1(",pc1_varlist3,"%)",sep=""))+ylab(paste("PC2(",pc2_varlist3,"%)",sep=""))+
  theme(axis.text=element_text(size=10),        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title=element_text(size=16,face="bold"))+
  geom_point(data=g1, colour="Black", shape=c(18), size=6) +
  theme(legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = 0)+
  theme(legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))+
  scale_x_continuous(limits = c(-45, 70))+
  scale_y_continuous(limits = c(-65, 50))+
  scale_colour_manual(name = "Collection site",values = group.col)+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1))
dev.off()

