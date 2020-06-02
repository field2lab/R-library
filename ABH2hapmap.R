ABH2hapmap <- function(file){
  hapmap <- read.csv(file, header=TRUE, row.names=1, sep=",",stringsAsFactors=FALSE)
  samples <- names(hapmap)[-1]
  loci <- row.names(hapmap)
  
  # set up conversion table
  s <- as.character(c("A","H","B","-"))
  
  ac<- c("A","M","C","N")
  ca= c("C","M","A","N")
  ag<- c("A","R","G","N")
  ga<- c("G","R","A","N")
  at<- c("A","W","T","N")
  ta<- c("T","W","A","N")
  cg<- c("C","S","G","N")
  gc<- c("G","S","C","N")
  ct<- c("C","Y","T","N")
  tc<- c("T","Y","C","N")
  gt<- c("G","K","T","N")
  tg<- c("G","K","T","N")
  
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
  names(conv) <- as.character(c("AC","CA","AG","GA","AT","TA","CG","GC",
                   "CT","TC","GT","TG"))
  
  # Pull out and convert genotypes
  S <- length(samples)
  SNPmatrix <- matrix(NA,length(loci),S)   # set up list of SNPbin objects
  
  total <- S
  # create progress bar
  pb <- txtProgressBar(min = 1, max = S, style = 3)
  
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=T, USE.NAMES=FALSE)
    mygen[sapply(mygen , is.null)] <- NA
    # create SNPbin object for this individual
    SNPmatrix[,i] <- unlist(mygen)
    setTxtProgressBar(pb, i)
  }
  finalhapmap=as.data.frame(SNPmatrix)
  colnames(finalhapmap)=samples
  finalhapmap[,c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode")] <- NA
  finalhapmap=finalhapmap[,c((length(finalhapmap)-10):length(finalhapmap),1:(length(finalhapmap)-11))]
  finalhapmap$alleles=hapmap$Type
  finalhapmap$alleles=paste(substr(finalhapmap$alleles,1,1),substr(finalhapmap$alleles,2,2),sep= "/")
  finalhapmap[,1]=loci
  return(finalhapmap)
}
mygeno=ABH2hapmap(file.choose())
write.table(mygeno,sep=",",eol = "\n",quote=FALSE,file = "SOAPhapmap.txt")
