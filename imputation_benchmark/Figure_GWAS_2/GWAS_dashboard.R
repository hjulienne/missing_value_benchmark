library(MASS)
library(ggplot2)
library(cowplot)
library(mice)
library(mclust)
library(MGMM)
library("missForest")
library('VIM')
library(FactoMineR)
library("factoextra")
options(warn=-1)

setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_GWAS_2/")
source("../Imputation_perf_functions.R")
source("../Adjusted_rand_index_table.R")

X = read.csv("./GWAS_cl.csv", row.names=1)
nsamp = dim(X)[1]

ncomponent = sum(grepl("z_", names(X)))
miss_rate = seq(0.0, 0.4, 0.05)
ncomponent
# Get cluster on complete data as reference
X["component"] = X[,"Cluster"]
Xtest  = insertNA(X, 0.2, 34)
opt_clust = 3
i=1
which.max(apply(X[,1:ncomponent],2,var,na.rm=TRUE))

for(i in 1:20)
{
  fo = paste0("./benchmark_data/Adj_RI_GWAS_",i,".csv")
  SR = as.numeric(substr(as.character(as.numeric(Sys.time())), 13,17))
  Adj_rand_index_table(X, ncomponent, opt_clust, miss_rate, fo=fo, seed=SR)#, fix.means=F,always.random=T)
}

#
# insertNA <- function(X, fraction){
#     Nna_bycol = ceiling(fraction * dim(X)[1])
#     for(j in 1:dim(X)[2])
#     {
#         ids = sample(dim(X)[1], Nna_bycol, replace=FALSE)
#         X[ids,j]=NA
#     }
#     return(X)
# }
#
# miss_rate = seq(0.0, 0.2, 0.01)
# ## Store results of different simulation
# Adjusted_Rand_index_table = data.frame( matrix(rep(0, 8*length(miss_rate)), ncol=8))
#
# row.names(Adjusted_Rand_index_table) = miss_rate
# names(Adjusted_Rand_index_table) = c("GMM_filtered", "GMM", "missForest", "KNN", "MICE","MICE_filtered", "mean", "median")




#
# for(mr in miss_rate){
#
#     Mt2=NULL
#
#     X_miss = X
#     X_miss[,1:ncomponent] = insertNA(as.matrix(X_miss[, 1:ncomponent]), mr)
#
#     Mt2 = compute_data_cluster_robust(as.data.frame(X_miss[,1:ncomponent]), 4)
#
#     X_miss["cluster"] = Mt2$clustered_data[,"cluster"]
#     X_miss["entropy"] = Mt2$clustered_data[,"entropy"]
#
#     if(!is.null(Mt2))
#     {
#       Adjusted_Rand_index_table[as.character(mr), "GMM"] = adjustedRandIndex(X_miss$component, X_miss$cluster)
#       # Rand index when we filter out ambiguous Assignation
#       well_identified = which(X_miss$entropy < 0.25)
#       Adjusted_Rand_index_table[as.character(mr), "GMM_filtered"] = adjustedRandIndex(X_miss[well_identified, ]$component, X_miss[well_identified, ]$cluster)
#     }
#     #######################################################
#     #######################################################
#     # Evaluation basic imputation strategies (mean, median)
#     #######################################################
#     X_mean = X_miss[,1:(ncomponent+1)]
#     X_median = X_miss[,1:(ncomponent+1)]
#
#     X_mean[,1:ncomponent] = impute_by_fun(X_miss[,1:ncomponent])
#     X_median[,1:ncomponent] = impute_by_fun(X_miss[,1:ncomponent], median)
#
#     Mt_mean <- compute_data_cluster_robust(X_mean[,1:ncomponent], 4)
#     Mt_median <- compute_data_cluster_robust(X_mean[,1:ncomponent], 4)
#     if(!is.null(Mt_mean)){Adjusted_Rand_index_table[as.character(mr), 'mean'] = adjustedRandIndex(Mt_mean$clustered_data[,"cluster"], X$component)}
#     if(!is.null(Mt_median)){Adjusted_Rand_index_table[as.character(mr), 'median'] = adjustedRandIndex(Mt_median$clustered_data[,"cluster"], X$component)}
#
#     #######################################################
#     #######################################################
#     # Evaluation of MICE imputation
#     # Let's trap the mice ;)
#     #######################################################
#     #######################################################
#
#     X_imputed = X_miss[,1:ncomponent]
#     row.names(X_imputed) = 1:nsamp
#     n_imp = 5
#     imputation = mice(X_imputed[,1:ncomponent], m=n_imp, maxit = 50, method = 'pmm', seed = 500)
#     # store cluster assignation for every imputation
#     cl_mice = data.frame(matrix(0,dim(X)[1], ncol=n_imp))
#     names(cl_mice) = 1:n_imp
#
#
#     for(i in 1:n_imp){
#         dat = complete(imputation, i)
#         Mt3 = compute_data_cluster_robust(dat, 4)
#         if(!is.null(Mt3))
#         {
#           sort_clust = rank(sapply(Mt3$MNMmix@Means, function(x){x[[1]]}))
#           cl_mice[,i] = sort_clust[Mt3$MNMmix@Assignments$Assignment]
#         }else{
#           print(paste0("Clustering failed for mice ",i))
#         cl_mice[i] = NULL
#         }
#     }
#
#     #X_imputed[,1:ncomponent] = complete(imputation)
#     X_imputed["cluster"] = apply(cl_mice, 1, getmode)
#     X_imputed["Proba_mode"] = apply(cl_mice[1:nsamp,]==X_imputed$cluster[1:nsamp], 1, mean)
#
#     Adjusted_Rand_index_table[as.character(mr), 'MICE'] = adjustedRandIndex(X_imputed$cluster, X$component)
#     Adjusted_Rand_index_table[as.character(mr), 'MICE_filtered'] = adjustedRandIndex(X_imputed[X_imputed$Proba_mode > 0.5,]$cluster, X[X_imputed$Proba_mode > 0.5,]$component)
#
#     ################################################################
#     # Evaluate missForest
#     ################################################################
#
#     RF_imp <- missForest(X_miss[,1:ncomponent])
#     X_RF = RF_imp$ximp[,1:ncomponent]
#     Mt4 <- compute_data_cluster_robust(X_RF[,1:ncomponent], 4)
#     if(!is.null(Mt4)){
#     X_RF["cluster"]= Mt4$MNMmix@Assignments$Assignment
#     #X_RF["entropy"] = apply(Mt4@Responsibilities[,-1], 1, shannon_entropy)
#     Adjusted_Rand_index_table[as.character(mr), 'missForest'] = adjustedRandIndex(X_RF$cluster, X$component)
#   }
#
#     ############################################################
#     # kNN
#     ############################################################
#     X_kNN <- kNN(X_miss[,1:ncomponent])
#     Mt5 <- compute_data_cluster_robust(X_kNN[, 1:ncomponent], 4)
#
#     X_kNN["cluster"]= Mt5$MNMmix@Assignments$Assignment
#     #X_kNN["entropy"] = apply(Mt5@Responsibilities[,-1], 1, shannon_entropy)
#     Adjusted_Rand_index_table[as.character(mr), 'KNN'] = adjustedRandIndex(X_kNN$cluster, X$component)
#
#     print(Adjusted_Rand_index_table)
#
#     write.csv(Adjusted_Rand_index_table, 'Adj_rand_index_GWAS.csv')
# }
# write.csv(Adjusted_Rand_index_table, 'Adj_rand_index_GWAS.csv')
