library(Matrix)
library(lme4)
library(MASS)
library(car)
library(reshape2)

## Read in Brix dataset
phenodata = read.csv(file.choose())

# fit the model
str(phenodata)
phenodata$Year=as.factor(phenodata$Year)
phenodata$Block=as.factor(phenodata$Block)


anova=aov(Transformed~ MixID+Block+Year+MixID:Year,data=phenodata)
summary(anova)
boxcox(anova,data=phenodata, lambda = seq(-10, 2, 1/10), 
       plotit = TRUE, eps = 1/50, xlab = expression(lambda),
       ylab = "log-Likelihood")

leveneTest(Transformed~MixID*Year, phenodata)
