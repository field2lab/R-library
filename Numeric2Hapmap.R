#####################################convert to genlight##############################
Numeric2hapmap <- function(file){
  require(adegenet)
  hapmapraw1 <- read.csv(file, header=TRUE,  sep="\t",
                         stringsAsFactors=FALSE,check.names=F)

  hapmap1=hapmapraw1[!(hapmapraw1$alleles=="A"),]
  hapmap1=hapmap1[!(hapmap1$alleles=="T"),]
  hapmap1=hapmap1[!(hapmap1$alleles=="G"),]
  hapmap1=hapmap1[!(hapmap1$alleles=="C"),]
  hapmap1[is.na(hapmap1)]="N"
  hapmap1[6:11]=NA
  
hapmap=hapmap1
hapmap[12:302]=apply(hapmap1[12:302], 2, as.character)
  
  samples <- names(hapmap[12:302])
  loci <- hapmap$alleles
  

  # set up conversion table
  s <- c("0","1","2","N")

  
  ac<- as.character(c("A","M","C","N"))
  ca= as.character(c("C","M","A","N"))
  ag<- as.character(c("A","R","G","N"))
  ga<- as.character(c("G","R","A","N"))
  at<- as.character(c("A","W","T","N"))
  ta<- as.character(c("T","W","A","N"))
  cg<- as.character(c("C","S","G","N"))
  gc<- as.character(c("G","S","C","N"))
  ct<- as.character(c("C","Y","T","N"))
  tc<- as.character(c("T","Y","C","N"))
  gt<- as.character(c("G","K","T","N"))
  tg<- as.character(c("G","K","T","N"))
  
  names(ac) <- s
  names(ag) <- s
  names(at) <- s
  names(cg) <- s
  names(ct) <- s
  names(gt) <- s
  names(ca) <- s
  names(ga) <- s
  names(ta) <- s
  names(gc) <- s
  names(tc) <- s
  names(tg) <- s
  conv <- list(ac,ca,ag,ga,at,ta,cg,gc,ct,tc,gt,tg)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G")
  
  S <- length(samples)
  SNPmatrix <- matrix(NA,length(loci),S)
  
  total <- S
  # create progress bar
  pb <- txtProgressBar(min = 1, max = S, style = 3)
  
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[2]], gen=hapmap[[11+i]],SIMPLIFY = T, USE.NAMES=FALSE)
    # create SNPbin object for this individual
    SNPmatrix[,i] <- mygen
    setTxtProgressBar(pb, i)
  }
  
  finalhapmap=as.data.frame(SNPmatrix,stringsAsFactors=FALSE)

  colnames(finalhapmap)=samples
  finalhapmap[,c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode")] <- NA
  finalhapmap=finalhapmap[,c((length(finalhapmap)-10):length(finalhapmap),1:(length(finalhapmap)-11))]
  finalhapmap$alleles=hapmap$alleles
  finalhapmap[,1]=hapmap$`rs#`
  finalhapmap[,3]=hapmap$chrom
  finalhapmap[,4]=round(hapmap$pos, digits = 0)
  finalhapmap[,5]="+"
  finalhapmap=finalhapmap[order(finalhapmap$chrom,finalhapmap$pos),]
return(finalhapmap)
}
mygeno=Numeric2hapmap(file.choose())
mygeno[,302]
write.table(mygeno,"hapmapforTassel.hmp.txt",row.names=F,sep="\t",quote = FALSE,eol = "\n")

myhapmap<- read.csv(file.choose(), header=TRUE,  sep="\t", stringsAsFactors=FALSE,check.names=F)
