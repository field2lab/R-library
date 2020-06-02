EmptyNA= read.csv(file.choose(),check.names = F)
EmptyNA[is.na(EmptyNA)]="."
SASfile=EmptyNA
write.csv(SASfile,file="SASfile.csv")
