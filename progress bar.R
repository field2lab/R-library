total <- 20
# create progress bar
pb <- txtProgressBar(min = 5, max = 21, style = 3)
for(i in 1:16){
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i+5)
}
close(pb)

Estimates=matrix(NA,255,1)
Estimates[,1]=temp[,1]
for (i in 5:21){
  inputData.entry =lmer(inputData[,i] ~ entry +(1|location)+(1|location:block)+ (1|entry:location), data=inputData)
  k=lsmeansLT(inputData.entry, test.effs="entry")
  temp=cbind(as.character(k[[1]]$entry),k[[1]]$Estimate,k[[1]][,3])
  Estimates=merge(Estimates, temp, by=1,all=TRUE)
  setTxtProgressBar(pb, i)
}