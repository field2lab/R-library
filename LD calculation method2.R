#Fitting method2###########
#452
n = 2*242
Cstart <- c(C=1)
dist15000=LDall15000$distance
rsq15000=LDall15000$r2
modelC15000 <- nls(r2 ~ ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))), data=LDall15000, start=Cstart, control=nls.control(maxiter=100))
rho15000 <- summary(modelC15000)$parameters[1]
newrsq15000 <- ((10+rho15000*LDall15000$distance)/((2+rho15000*LDall15000$distance)*(11+rho15000*LDall15000$distance)))*(1+((3+rho15000*LDall15000$distance)*(12+12*rho15000*LDall15000$distance+(rho15000*LDall15000$distance)^2))/(n*(2+rho15000*LDall15000$distance)*(11+rho15000*LDall15000$distance)))
newfile15000 <- data.frame(LDall15000$distance, newrsq15000)
maxld15000 <- max(newfile15000$newrsq15000)
halfdecay15000 = maxld15000*0.5
halfdecaydist15000<- newfile15000$LDall15000.dist[which.min(abs(newfile15000$newrsq15000-halfdecay15000))]
newfile15000 <- newfile15000[order(newfile15000$LDall15000.dist),]
f115000 <- data.frame(newfile15000$LDall15000.dist, newfile15000$newrsq15000)
xval15000 <- f115000[which.min(abs(0.2 - f115000$newfile15000.newrsq15000)),] #find x value where y=0.2
xval15000[1,1]


################## create the plot###########################
#png(file="LDall70000.png",width=12.57703333,height=8,units="in",res=700)
gaps=seq(from=1, to=82118,by=500)
line1dist=newfile15000$LDall15000.dist[gaps[c(1,80,110,132,142:length(gaps))]]
line1rsq=newfile15000$newrsq15000[gaps[c(1,80,110,132,142:length(gaps))]]
plot(1, type="n", xlab="distance (bp)", ylab="LD (r^2)", xlim=c(0, 5000), ylim=c(0, 0.6))
#(newfile40000$LDall40000.dist, cex = .1,pch=0, newfile40000$newrsq40000, col="black")
#points(newfile30000$LDall30000.dist, cex = .1, pch=1, newfile30000$newrsq30000, col="black")
#points(newfile20000$LDall20000.dist, cex = .1, pch=2, newfile20000$newrsq20000, col="black")
lines(line1dist, type="o",line1rsq, col="black")
#points(newfile10000$LDall10000.dist, cex = .1,pch=4,  newfile10000$newrsq10000, col="black")
lines(frame10000, col="blue")
lines(frame15000, col="red")
lines(frame20000, col="orange")
lines(frame30000, col="green")
lines(frame40000, col="purple")
abline(h=0.2, col="black")
legend(3000, 0.6, legend=c("10000 with method1", "15000 with method1","20000 with method1","30000 with method1","40000 with method1"),
       lty=1,col=c("blue", "red", "orange", "green", "purple"),
       title="Distance(bp)",cex=0.8)
legend(3000, 0.3, legend=c("15000 with method2"),
       lty=1,pch=1,col=c("black"),cex=0.8)
abline(v=xval[1,1], col="green")
mtext(round(xval[1,1],2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
#dev.off()