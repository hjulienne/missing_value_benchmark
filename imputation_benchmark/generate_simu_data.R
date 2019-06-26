library(MASS)
library(ggplot2)
library(cowplot)
library(mice)
library(mclust)
library(MGMM)
library("missForest")
library('VIM')
library("colourlovers")
setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/")
#source('../clustering_analysis_function.R')
source("./Imputation_perf_functions.R")

###
# Create a bivariate dataset from a gaussian mixture with 4 components
# The component are placed so knowledge of the first variables should be enough
# to retrieve while the second variable alone shouldn't be enough
# The Covariance matrix is chosen so the local correlation inside one component
# goes in a different direction than the global correlation
###

Sigma = list()
means = list()
for(i in 1:4){
  means[[i]] = sample(-5:5, 2)
  covxy = runif(1, -0.9, 0.9)
  Sigma[[i]] = matrix(c(0.9, covxy, covxy, 0.9),2,2)
}
=
X = rGMM(2000, d=2, k=4, M = means, S=Sigma)
component = row.names(X)
X = as.data.frame(X)
names(X) = c("y1", 'y2')
X["component"] = component

write.csv(X, "data_bivariate_simulation.csv")
