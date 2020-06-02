ACBD_generator=function(ncheck,ntrt,nperB,nrowperB,outputdir){
#Enviroment and parameter setup
  setwd(outputdir)
  is.even <- function(x) x %% 2 == 0
  #is.odd <- function(x) x %% 2 != 0
  nB=ceiling(ntrt/(nperB-ncheck))
  nperrow=nperB/nrowperB
  ntrtperB=nperB-ncheck
  checkvector1=as.data.frame(sample(1:ncheck,replace=FALSE))
  checkvector2=as.vector(paste("C",checkvector1[,1],sep=""))
  temp1=c(sample(1:ntrt, ntrt, replace=F))
  temp1new=c(temp1,rep(NA, length(1:(nB*nperB-(ntrt+nB*ncheck)))))
  temp4=as.data.frame(matrix(NA,nrowperB*nB,nperrow))
#ACBD entry map
for(i in 1:(nB-1)){
temp2=as.data.frame(as.character(sample(c(temp1new[((i-1)*ntrtperB+1):(i*ntrtperB)],checkvector2),replace=F)))
colnames(temp2)="X"
temp3=as.data.frame(matrix(NA,nrowperB,nperrow))
for (j in 1:nrowperB){
temp3=rbind(as.character(temp2[((j-1)*nperrow+1):(j*nperrow),1]),temp3)  
}
temp4=rbind(temp3[1:nrowperB,],temp4)
}

temp3last=as.data.frame(matrix(NA,nrowperB,nperrow))
temp2last=as.data.frame(c(rep(NA, length(1:(nB*nperB-(ntrt+nB*ncheck)))),
                          as.character(sample(c(temp1new[((nB-1)*ntrtperB+1):ntrt],checkvector2),replace=F))))
for (j in 1:nrowperB){
  if (is.even(j)) temp3lastx=rev(as.character(temp2last[((j-1)*nperrow+1):(j*nperrow),1])) else temp3lastx=as.character(temp2last[((j-1)*nperrow+1):(j*nperrow),1]) 
  temp3last=rbind(temp3last,temp3lastx)
}
if (is.even(nB*nrowperB-1)) temp3lastnew=rev(temp3last[(nrow(temp3last)-nrowperB+1):nrow(temp3last),]) else temp3lastnew=temp3last[(nrow(temp3last)-nrowperB+1):nrow(temp3last),] 
temp4=rbind(temp3lastnew,temp4)

tempmap1=temp4[1:(nrowperB*nB),]
tempmap1[is.na(tempmap1)]="filler"
map1=as.matrix(tempmap1)
rownames(map1)=NULL
rownames(map1)[1:nrow(map1)]=c(' ')
rownames(map1)[seq(nrowperB,nB*nrowperB,by=nrowperB)]=paste("Block",c(nB:1))
colnames(map1)=paste("Path",c(1:nperrow))

#ACBD field map
Plotseq=as.data.frame(matrix(NA,nrowperB*nB,nperrow))
for(i in 1:(nB*nrowperB)){
  #Plotseq1=t(as.data.frame(seq(((i-1)*nperrow+101*i),i*nperrow+100*i,by=1)))
  Plotseq1=t(as.data.frame(seq((i*100)+1,(i*100)+nperrow,by=1)))
  Plotseq[i,]=if(is.even(i)) apply(Plotseq1,1,rev) else Plotseq1
}
Plotseq=as.matrix(Plotseq)
map2=apply(Plotseq,2,rev)
rownames(map2)=NULL
rownames(map2)[1:nrow(map1)]=c(' ')
rownames(map2)[seq(nrowperB,nB*nrowperB,by=nrowperB)]=paste("Block",c(nB:1))
colnames(map2)=paste("Path",c(1:nperrow))

#Entry list
temp5combine=as.data.frame(matrix(NA,1,0))
rownames(tempmap1)=NULL
for(i in 1:nrow(map1)){
temp5=tempmap1[nB*nrowperB-i+1,]
if (is.even(i)) temp5new=rev(temp5) else temp5new=temp5
temp5combine=cbind(temp5combine,temp5new)
}

tempentrylist=t(as.matrix(temp5combine))
entrylist=as.data.frame(tempentrylist)
spots=seq(1,length(tempentrylist),by=nperB)
for(i in 1:length(spots)){
  entrylist[spots[i]:(spots[i]+nperB-1),2]=i
}
rownames(entrylist)=NULL
#entrylist[,3]=seq(101,100+nB*nperB,by=1)
for(i in 1:(nB*nrowperB)){
entrylist[((i-1)*nperrow+1):(i*nperrow),3]=seq((i*100)+1,(i*100)+nperrow,by=1)
}
colnames(entrylist)=c("Entry","Block","Plot")
entrylist=entrylist[,c("Plot","Block","Entry")]

#combined map

map2new=as.matrix(map2)
map1new=as.matrix(map1)

combinedmap=matrix(paste(map2new, map1new, sep= "_"), 
                    nrow=nrow(map2new))
colnames(combinedmap)=colnames(map2)
rownames(combinedmap)=rownames(map2)

#write csv files
data_list=list(map1, map2,combinedmap,entrylist)
names(data_list)=c("Entry map","Field map","combined map","Entry list")

return(for(i in names(data_list)){write.csv(data_list[[i]], paste0(i,".csv"))})
}

ACBD_generator(3,147,24,1,"/Users/Rorshach/Desktop")
ACBD_generator(3,153,25,1,"C:/Users/Jia/Desktop")
ACBD_generator(3,193,24,1,"C:/Users/Jia/Desktop")
ACBD_generator(3,193,27,9,"C:/Users/Jia/Desktop")

