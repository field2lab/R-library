library(sp)
library(rgdal)
library(PBSmapping)
library(maptools)
library(raster)
library(rgeos)
library(plyr)
library(ggplot2)
library(reshape2)
library(maps)
library(mapproj)
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}
################################################################################
# Read ESRI Shape File into SpatialPolygonsDataFrame
# rgdal readOGR() method gets projection information.
Hardiness <- readOGR("/home/rorshach/Desktop/Hardinessmap/ophz-c", "ophz-c")
Hardiness2<- spTransform(Hardiness, CRS("+proj=longlat +datum=WGS84"))
################################################################################
#EcoLegend=colsplit(EcoRegion2@data$L3_KEY," ",c("Code","name"))
#EcoLegend2=ddply(EcoLegend, ~ name, summarise, Code= min(Code))
################################################################################
#EcoRegion2@data$id = rownames(EcoRegion2@data)
#EcoRegion2.points = fortify(EcoRegion2, region="id")
#EcoRegion2.df = join(EcoRegion2.points, EcoRegion2@data, by="id")
###############################################################################
HDcoords=as.data.frame(cbind(newDAPClist[,colnames(newDAPClist)=="Longitude"],newDAPClist[,colnames(newDAPClist)=="Latitude"]))
for (i in 1:nrow(HDcoords)){
  if(is.na(HDcoords[i,1])) newDAPClist$Hardinessz[i]="Unknown"
  else newDAPClist$Hardinessz[i]=as.character(over(SpatialPoints(HDcoords[i,],proj4string=CRS(proj4string(Hardiness2))), Hardiness2)$zone)
  print(i)
}
HDcoords2=as.data.frame(cbind(newDAPClist[,colnames(newDAPClist)=="Longitude"],newDAPClist[,colnames(newDAPClist)=="Latitude"]))
for (i in 1:nrow(HDcoords2)){
  if(is.na(HDcoords2[i,1])) newDAPClist$Hardinessz[i]="Unknown"
  else newDAPClist$HDtemp[i]=as.character(over(SpatialPoints(HDcoords2[i,],proj4string=CRS(proj4string(Hardiness2))), Hardiness2)$temp)
  print(i)
}
# outputPOPlist
###############################################################################
Poplistfixed=newDAPClist[!duplicated(newDAPClist$State2),]
Poplistfixed$Latitude1=as.character(dd2dms(Poplistfixed$Latitude,NS=T))
Poplistfixed$Longitude1=as.character(dd2dms(Poplistfixed$Longitude,NS=F))
write.csv(Poplistfixed,file = "97Poplist.csv" )


