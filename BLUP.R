library(Matrix)
library(lme4)
library(MASS)

## Read in Brix dataset
phenodata = read.csv(file.choose())

## Examine distribution of brix data
hist(Brix, col="gold")
boxplot(Brix~Loc, xlab="Location", ylab="Degrees Brix", main="Degrees Brix by Location", col="pink")

## BLUPS
# fit the model
str(phenodata)
anova1=aov(HD~ family/ID2+family+Year,data=phenodata)
boxcox(anova1,data=phenodata, lambda = seq(-10, 2, 1/10), 
       plotit = TRUE, eps = 1/50, xlab = expression(lambda),
       ylab = "log-Likelihood")
phenodata$newgplant=sqrt(phenodata$Plant.Wt.)
phenodata$newTN.no.plant=sqrt(phenodata$Tiller.No.)

phenomodel1= lmer(Plant.Wt.~ 1+(1|SCA)+(1|Year),data=phenodata)
phenomodel2= lmer(Tiller.Wt.~ 1+(1|Code)+(1|Year),data=phenodata)
phenomodel3= lmer(Stem.D.~ 1+(1|Code)+(1|Year),data=phenodata)
phenomodel4= lmer(newTN.no.plant~ 1+(1|Code)+(1|Year),data=phenodata)
phenomodel5= lmer(Tiller.L.~ 1+(1|Code)+(1|Year),data=phenodata)
phenomodel6= lmer(Leaf.No.~ 1+(1|Code)+(1|Year),data=phenodata)

BLUPsall=matrix(,nrow=248)
# estimate BLUPs$Repeatabilities
g.plantBLUP = ranef(phenomodel1,,condVar=T)
g.plantCodeBLUP = g.plantBLUP$SCA
dotplot(g.plantBLUP, scales = list(x = list(relation = 'free')))[["SCA"]]

fit1= lmer(newgplant~ (1|Code)+(1|family)+Year,data=phenodata)
vc1 = VarCorr(fit1)
Genetic_var1 = attr(vc1$Code,'stddev')[1]^2+attr(vc1$family,'stddev')[1]^2
Other_var1=attr(vc1,'sc')^2
R1= Genetic_var1 /(Genetic_var1+Other_var1)

g.tillerBLUP = ranef(phenomodel2)
g.tillerCodeBLUP = g.tillerBLUP$Code
fit2= lmer(Tiller.Wt.~ (1|Code)+(1|family)+Year,data=phenodata)
vc2 = VarCorr(fit2)
Genetic_var2 = attr(vc2$Code,'stddev')[1]^2+attr(vc2$family,'stddev')[1]^2
Other_var2=attr(vc2,'sc')^2
R2= Genetic_var2 /(Genetic_var2+Other_var2)

D.averageBLUP = ranef(phenomodel3)
D.averageCodeBLUP = D.averageBLUP$Code
fit3= lmer(Stem.D.~ (1|Code)+(1|family)+Year,data=phenodata)
vc3 = VarCorr(fit3)
Genetic_var3 = attr(vc3$Code,'stddev')[1]^2+attr(vc3$family,'stddev')[1]^2
Other_var3=attr(vc3,'sc')^2
R3= Genetic_var3 /(Genetic_var3+Other_var3)

TN.no.plantBLUP = ranef(phenomodel4)
TN.no.plantCodeBLUP = TN.no.plantBLUP$Code
fit4= lmer(newTN.no.plant~ (1|Code)+(1|family)+Year+(1|family:Year),data=phenodata)
vc4 = VarCorr(fit4)
Genetic_var4 = attr(vc4$Code,'stddev')[1]^2+attr(vc4$family,'stddev')[1]^2
Other_var4=attr(vc4,'sc')^2
R4= Genetic_var4 /(Genetic_var4+Other_var4)

PH.mBLUP = ranef(phenomodel5)
PH.mCodeBLUP =PH.mBLUP$Code
fit5= lmer(Tiller.L.~ (1|Code)+(1|family)+Year+(1|family:Year),data=phenodata)
vc5 = VarCorr(fit5)
Genetic_var5 = attr(vc5$Code,'stddev')[1]^2+attr(vc5$family,'stddev')[1]^2
Other_var5=attr(vc5,'sc')^2
R5= Genetic_var5/(Genetic_var5+Other_var5)

LN.no.tillerBLUP = ranef(phenomodel6)
LN.no.tillerCodeBLUP =LN.no.tillerBLUP$Code
fit6= lmer(Leaf.No.~ (1|Code)+(1|family)+Year+(1|family:Year),data=phenodata)
vc6= VarCorr(fit6)
Genetic_var6 = attr(vc6$Code,'stddev')[1]^2+attr(vc6$family,'stddev')[1]^2
Other_var6=attr(vc6,'sc')^2
R6= Genetic_var6/(Genetic_var6+Other_var6)

BLUPsall=cbind(g.plantCodeBLUP,g.tillerCodeBLUP,D.averageCodeBLUP,TN.no.plantCodeBLUP,PH.mCodeBLUP,LN.no.tillerCodeBLUP)
colnames(BLUPsall)=c("Plant Wt.","Tiller Wt.","Stem D.","Tiller No.","Tiller L.","Leaf No.")
Repit=rbind(R1,R2,R3,R4,R5,R6)
rownames(Repit)=c("Plant Wt.","Tiller Wt.","Stem D.","Tiller No.","Tiller L.","Leaf No.")
colnames(Repit)=c("Repeatability")




# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="BrixLineBLUPS.csv")

## Creating plots with the BLUPs
# Create a numeric vector with the BLUP for each line
IDBLUP = g.plantIDBLUP[,1]
# Create a histogram with the BLUP for each line
hist(IDBLUP, col="brown")

## Compare BLUP to line averages on a scatterplot
lmean = tapply(phenodata$g.plant, phenodata$ID, na.rm=T, mean)
plot(IDBLUP, lmean)

