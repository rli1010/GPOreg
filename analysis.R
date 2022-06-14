library(Rcpp) 
library(truncnorm)

#simda=read.csv("simdata.csv")
#source("functions/utils.R")

trans.reg(cbind(Y1,Y2,Y3,Y4) ~ Z1+Z2+Z3+Z4, data = simda)
