library(Rcpp) 
library(truncnorm)

simda<- read.csv("simdata.csv",header=TRUE)
sourceCpp("code.cpp")
source("functions.R")

##Fitting start here
##Data rows with missing data need to be removed before fitting

#Transformed linear model
ff=trans.reg(cbind(Y1,Y2,Y3,Y4) ~ Z1+Z2+Z3+Z4, data = simda)
ff

#Standardized to unit norm
ff$Coef[-1]/sqrt(sum(ff$Coef[-1]^2))


#MRC
#Set resamp to an integer>=200 for real application

#As the sign of beta[1] is fixed to be positive, set the first variable in varlist as a covariate with positive coef (e.g.,  a significant variable with positive 
#coef in the transformed linear model)

mono.fit=mono.reg(Ylist=c("Y1","Y2","Y3","Y4"),varlist=c("Z1","Z2","Z3","Z4"),
                  da=simda,resamp=30)  

mono.fit$est
