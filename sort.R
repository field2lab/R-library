sortdata=read.csv(file.choose())
workdata=as.data.frame(sortdata)
workdata$T3=as.numeric(workdata$T3)
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
newdata=Nth.delete(workdata, 3)
rownames(newdata)=NULL
#workdata[7:16]<- as.numeric(as.matrix(workdata[7:16]))

newsort2=NULL
colnames(newsort2)=colnames(cbind(newdata[1,c(1:15)],newdata[2,7:16],newdata[3,7:16]))
for (i in 1:(nrow(newdata)/2)){
  newsort1=cbind(newdata[3*i-2,c(1:16)],newdata[3*i-1,7:16],newdata[3*i,7:16])
  newsort2=rbind(newsort2,as.data.frame(newsort1[1,]))
}
sorteddata=newsort2[-1,-6]
rownames(sorteddata)=NULL
colnames(sorteddata)=c("Person","Name",	"Plot",	"Location",	"Date",	"I1R1",	"I1R2",	"I1R3",	"I1R4",	"I1R5",	"I1R6",	"I1R7",	"I1R8",	"I1R9",	"I1R10","I2R1",	"I2R2",	"I2R3",	"I2R4",	"I2R5",	"I2R6",	"I2R7",	"I2R8",	"I2R9",	"I2R10","I3R1",	"I3R2",	"I3R3",	"I3R4",	"I3R5",	"I3R6",	"I3R7",	"I3R8",	"I3R9",	"I3R10")

write.csv(sorteddata, file="sorteddata.csv")
