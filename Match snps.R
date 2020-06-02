Mappedmp=read.csv(file.choose(),header=T,sep="")
library(stringr)
Mappedmp[,1]=sapply(str_split(Mappedmp[,1],"_"),'[',1)
Pmap=read.csv(file.choose(),header=T)

toDelete=seq(0, nrow(Pmap), 2)
Pmap=Pmap[-toDelete, ]
Pmap=Pmap[!is.na(Pmap[,2]),]

Tasslefile=Mappedmp[Mappedmp[,1]%in%Pmap[,1],]
colnames(Tasslefile)[1]="rs"

Tasslefile[,3]=Pmap[match(Tasslefile[,1],Pmap[,1]),2]
Tasslefile[,4]=Pmap[match(Tasslefile[,1],Pmap[,1]),3]
rownames(Tasslefile)=NULL
write(Tasslefile, file="Tasselehmp.txt", sep = "\t")
write.table(Tasslefile,file = "Tasselehmp.txt",row.names=F,quote=F, sep = " ")
str(Tasslefile)
