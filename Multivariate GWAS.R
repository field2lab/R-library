library(sommer)
####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
####=========================================####
####=========================================####
#### EXAMPLE 1
#### GWAS in diploids
####=========================================####
####=========================================####
data(CPdata)
head(CPpheno)
CPgeno[1:4,1:4]
#### create the variance-covariance matrix
A <- A.mat(CPgeno)
#### look at the data and fit the model
head(CPpheno)
mix1 <- GWAS2(cbind(color,Yield)~1,
random=~g(id),
rcov=~units,
G=list(id=A),
M=CPgeno,
data=CPpheno)
summary(mix1)
