myfq=read.csv(file.choose())
newfq=as.vector(gsub("[][/]","",myfq$sequence))
myfq$newseq=newfq

myname=paste(">",myfq$name,sep="")
myfq$newname=myname
myfq=myfq[myfq$newseq!="Unavailable",]

t2=c()
for (i in 1: nrow(myfq)){
  t1=t(cbind(myfq[i,]$newname,myfq[i,]$newseq))
  t2=rbind(t1,t2)
}

myfq$name=as.character(myfq$name)
library(seqinr)
write.fasta(as.list(myfq$newseq),c(myfq$name),nbchar = 10000, file.out = "consensus_oat.fasta")

myread=read.fasta(file.choose(),as.string = T)
