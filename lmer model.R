library(Matrix)
library(lme4)
library(MASS)
library(lmerTest)

## Read in Brix dataset
phenodata = read.csv(file.choose())

fit1= lmer(DM_Mg_ha~ Variety.No.+(1|Block..Rep..No.),data=phenodata)
anova(fit1)



phenodata$Variety.No.=as.factor(phenodata$Variety.No.)
phenodata$Block..Rep..No.=as.factor(phenodata$Block..Rep..No.)

summary(phenodata$Variety.No)


lsmeans(fit1)

lsmeans(fit1, pairwise~"Variety.No.", adjust="tukey")
