library(elevatr)
#cat("mapzen_key=mapzen-XXXXXXX\n", file = file.path(normalizePath("~/"), ".Renviron"), 
 #   append = TRUE)

set.seed(65.7)
examp_df <- data.frame(x = runif(10, min = -73, max = -71), y = runif(10, min = 41, 
                                                                      max = 45))
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Create and example data.frame with additional columns
cats <- data.frame(category = c("H", "H", "L", "L", "L", "M", "H", "L", "M", 
                                "M"))

examp_df2 <- data.frame(examp_df, cats)

# Create an example SpatialPoints
examp_sp <- SpatialPoints(examp_df, proj4string = CRS(prj_dd))

# Create an example SpatialPointsDataFrame
examp_spdf <- SpatialPointsDataFrame(examp_sp, proj4string = CRS(prj_dd), data = cats)

df_elev_epqs <- get_elev_point(examp_df, prj = prj_dd, src = "epqs")

data.frame(df_elev_epqs)


myelevdata=as.data.frame(cbind(newDAPClist[,7],newDAPClist[,6]))
myelev <- get_elev_point(myelevdata, prj = prj_dd, src = "epqs")
myelevdata=data.frame(myelev)


Qvaldata=cbind(myPCAdata,myelevdata$elevation, fastcladenew[,1:3])
colnames(Qvaldata)[16]="elevation"
Qvaldata$Ploidy=as.factor(Qvaldata$Ploidy)
Qvaldata$deme=as.factor(Qvaldata$deme)
Qvaldatanew=Qvaldata[Qvaldata$State2!=c("ME1","NE4","NE5"),]
Qvaldatanew=Qvaldatanew[Qvaldatanew$State2!=c("ND"),]
rownames(Qvaldatanew)=NULL
Qvaldatanew[,8]=as.character(Qvaldatanew[,8])
Qvaldatanew[,9]=as.character(Qvaldatanew[,9])
Qvaldatanew[92,8]='Red River'
Qvaldatanew[92,9]='Red River'
gallQ <- subset(Qvaldatanew, State == "Red River")
gall_1Q <- subset(Qvaldatanew, State2 == "NY2")
gall_2Q <- subset(Qvaldatanew, State2 == "NY1")
gall_1Q[1,9]='KST'
gall_2Q[1,9]='STP'
gallnewQ=rbind(gall_1Q,gall_2Q,gallQ)
gallnewQ$State2=factor(gallnewQ$State2, levels=c('KST','STP','Red River'), ordered = TRUE)
Qvaldatanew=Qvaldatanew[-c(1,2,92),]


lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(r)^2~"="~r2*","~~italic(P)~"<"~pv, 
                   #italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"<"~pv,
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        pv=0.001))
  as.character(as.expression(eq));                 
}

regdata=as.data.frame(cbind(Qvaldatanew$elevation,Qvaldatanew$fast1))
colnames(regdata)=c("x","y")


png(file="ECOdata2.png",width=12.57703333,height=8,units="in",res=700)
q=ggplot(Qvaldatanew,aes(x=elevation, y=fast1))+
  geom_jitter(size=5,stroke = 0.5, aes(shape=as.factor(Ploidy)))+theme_bw()+xlab("Elevation(m)")+ylab(paste("Posterior probability"))+
  #geom_text(aes(label=State2),hjust=0, vjust=-1.2)+
  geom_text(size=8,x = 200, y = 1.2, label = lm_eqn(regdata), parse = TRUE)+
  geom_smooth(method='lm',se=F,weight=1,color="black")+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(1,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())+
  
  guides(shape = guide_legend(order=2))+
  geom_point(data=gallnewQ,aes(x=gallnewQ$elevation,y=gallnewQ$fast1,color=State2),fill=c("orange","black","black"),shape=c(7,9,18),size=c(6,5,8),stroke = 1.5) +
  #scale_fill_manual(name = "Demes",values = c("Blue","Red"), labels=c("Deme1", "Deme2"))+
  #guides(fill = guide_legend(override.aes = list(shape=18,color =c("Red","Blue"),size=8),order=1))+
  scale_shape_manual(name = "Ploidy level",values = c(22,24,21), labels=c("2n=4x","2n=6x","2n=8x"))+
  guides(shape = guide_legend(override.aes = list(size=5,stroke=1.5),order=2))
q=q+scale_color_manual(name="Released varieties",values =c("black","black","orange"))+
  guides(color = guide_legend(override.aes = list(shape=c(7,9,18),size=c(6,5,8),color =c("black","black","orange"))),order=3)
q 
dev.off()