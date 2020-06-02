RCB_generator=function(nrep,nrow,entries,Plot_name,Entries_category,output_address){
  block=as.data.frame(matrix(NA,nrow,2*nrep))
  x=ifelse(nrow<=99,2,3)
  for(i in 1:nrep){
    block[,(2*i-1)]=sample(((i*10^x)+1):((i*10^x)+nrow),replace=FALSE) 
    colnames(block)[2*i-1]=Plot_name
  }
  for(j in 1:nrep){
    block[,2*j]=sample(entries,replace=FALSE)
    colnames(block)[2*j]=Entries_category
  }
  Field_design=block[sort.list(block[,1]),]
  Field_design[,2]=entries
  return(write.csv(Field_design,row.names =FALSE,file=output_address))
}

RCB_generator(4,12,c("IL102","PCG109","OC","SIDENY","SOYBEAN","M*G","KANLOW","CIR","20-102","IL99","IL104","17-116"),"Plot","Pops","C:/Users/Jia/Desktop/RCBDnursery.csv")



Split_plot=function(nrep,treatment_number,entries,treatment_name,output_address){
  x=ifelse(treatment_number*length(entries)<=99,2,3)
  block=as.data.frame(matrix(NA,treatment_number*length(entries),3*nrep))
  for(i in 1:nrep){
    m=sample(1:treatment_number,treatment_number, replace=FALSE)
    block[,(3*i-2)]=(((i*10^x)+1):((i*10^x)+treatment_number*length(entries)))
    colnames(block)[3*i-2]="Plot"
    colnames(block)[3*i-1]="Treatment"
    for (k in 1 :treatment_number){
      n=m[k]
      block[(length(entries)*k-(length(entries)-1)):(k*length(entries)),(3*i-1)]=paste(treatment_name,n,sep=" ")
    }}
  for(j in 1:nrep){
    colnames(block)[3*j]="Pops"
    for (y in 1:treatment_number) {
      block[(length(entries)*y-(length(entries)-1)):(y*length(entries)),3*j]=sample(entries,replace=FALSE)
    }}
  return(write.csv(block,row.names =FALSE,file=output_address))
}

Split_plot(4,4,c("IL102","PCG109","OC","M*G"),"HT","C:/Users/Jia/Desktop/SplitBnursery.csv")
