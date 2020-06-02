#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("ShortRead", version = "3.8")

myfq<- read.table("D:/Dropbox (UFL)/wheat/Dipendra/sequences/markerSeqs.txt",sep="\t")
DipG <- read.csv("D:/Dropbox (UFL)/wheat/Dipendra/2018/3_8_19/hete_removed_noindel_50missing_imputed.hmp.txt", header=T, sep="\t",
                 stringsAsFactors=FALSE)
is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0
myfqdata=myfq[is.even(as.numeric(rownames(myfq))),]
myseq=as.data.frame(myfqdata)

myfqdata2=myfq[is.odd(as.numeric(rownames(myfq))),]
mysnp=as.data.frame(myfqdata2)

selectedsnp=read.csv(file.choose(),sep=",")

mysnp[,1]=gsub(">chr","S",mysnp[,1])

selectedsnp$allele=DipG1[match(selectedsnp$SNP,DipG1$V1),2]
selectedsnp$seq=myseq[match(selectedsnp$SNP,mysnp[,1]),1]
selectedsnp$seqmark=paste(substr(selectedsnp$seq, 1, 199), "[",substr(selectedsnp$seq,200,200),"]" ,substr(selectedsnp$seq, 201, 400), sep = "")

write.csv(selectedsnp,"selectedsnp_with_seqs_4_23_19.csv",row.names = F)
