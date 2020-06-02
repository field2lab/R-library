library(LDcorSV)
data(data.test)
as.data.frame(data.test)
Geno<-data.test[[1]]
Geno.test<-Geno[,1:20]
V.WAIS<-data.test[[2]]
S.2POP<-data.test[[3]]

LD<-LD.Measures(Geno.test,V=V.WAIS,S=S.2POP,data ="G",supinfo=TRUE,na.presence=TRUE)
LD
plot(LD$r2,LD$r2vs)
plot(LD$r2,LD$r2v)
plot(LD$r2,LD$r2s)

dim()