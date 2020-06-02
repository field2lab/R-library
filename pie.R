library(ggplot2)
library(maps)
library(plotrix)
##########################################################################################################################################################
#plot map
all_states <- map_data("state")
states <- subset(all_states, region %in% c( "illinois",
                                            "indiana",
                                            "iowa","minnesota",
                                            "missouri", "north dakota",
                                            "south dakota", "wisconsin",
                                            "kansas",
                                            "oklahoma",
                                            "new york",
                                            "connecticut",
                                            "new hampshire",
                                            "massachusetts",
                                            "maine",
                                            "louisiana",
                                            "nebraska",
                                            "arkansas",
                                            "ohio",
                                            "pennsylvania"
                                           ))
p <- ggplot()
p <- p + geom_polygon( data=states, aes(x=long, y=lat, group = group),colour="black", fill="white" )
p



gps=mylist #Read input files
names(PPtable)[names(PPtable)=="1"]="c1"
names(PPtable)[names(PPtable)=="2"]="c2"
names(PPtable)[names(PPtable)=="3"]="c3"
names(PPtable)[names(PPtable)=="4"]="c4"
admix=as.data.frame(PPtable)

map("state") #Plot maps
map.axes() #Add axes
points(gps$Lon, gps$Lat,
       cex = admix$Num/5, col='red', pch=19) #To add just points
for (x in 2:nrow(gps)) { #To add arrows for migration models
  arrows(gps$Lon[1],gps$Lat[1],gps$Lon[x],
         gps$Lat[x],lty="dashed",code=2)
} #To add admixture plots - here I used K = 2.


map("state") 
map.axes()
for (x in 1:nrow(gps)) {
  floating.pie(gps$Longitude[x],gps$Latitude[x],c(admix$c1[x],admix$c2[x],
              admix$c3[x],admix$c4[x]),radius=0.4,
              col=c("red","blue","orange","black"))
  print(x)
  }
floating.pie(gps$Longitude[1],gps$Latitude[1],
             c(admix$c1[1],admix$c2[1],admix$c3[1],admix$c4[1]),radius=1,
              col=c("red","blue","orange","black"))



