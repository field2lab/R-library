#####################################convert to genlight##############################
ABH2genlight2 <- function(file){
  require(adegenet)
  ABH <- read.csv(file, header=TRUE, row.names=1, sep=",",
                         stringsAsFactors=FALSE)
  samples <- names(ABH)[-1]
  loci <- row.names(ABH)
    map = setNames(c(0,1,2,NA), c("A","H","B","-"))

  S <- length(samples)
  SBlist <- vector(mode="list",S)
  total <- S
  pb <- txtProgressBar(min = 1, max = S, style = 3)
  for(i in 1:S){
    mygen <- mapply(function(gen) map[unlist(gen)],
                    gen=ABH[[i+1]],
                    SIMPLIFY=T, USE.NAMES=FALSE)
    mygen[sapply(mygen , is.null)] <- NA
    SBlist[[i]] <- new("SNPbin", unlist(mygen))
    setTxtProgressBar(pb, i)
  }
  x <- new("genlight", SBlist,n.loc=length(loci),ind.names=samples,loc.names=loci)
  output=list(x)
  return(output)
}
mygeno=ABH2genlight2(file.choose())

hapmap=read.csv(file.choose(), header=TRUE,check.names=F, row.names=1, sep=",",
                stringsAsFactors=FALSE)
SOAPhapmap=hapmap
SOAPhapmap[,c("alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode")] <- NA
SOAPhapmap=SOAPhapmap[,c(290:299,1:289)]
SOAPhapmap$alleles=SOAPhapmap$Type
SOAPhapmap=SOAPhapmap[,-11]
SOAPhapmap$alleles=paste(substr(SOAPhapmap$alleles,1,1),substr(SOAPhapmap$alleles,2,2),sep= "/")

write.table(SOAPhapmap,sep="\t",file = "SOAPhapmap.txt")
SOAPsnp=as.data.frame(mygeno)