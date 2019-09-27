# A big function to compute the MGMM performances
# On data with missing values versus imputed data
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
source("~/Project/missing_value_benchmark/imputation_benchmark/clustering_analysis_function.R")

insertNA <- function(X, fraction, seed){
    Nna_bycol = ceiling(fraction * dim(X)[1])
    print(Nna_bycol)
    set.seed(seed)
    for(j in 1:dim(X)[2])
    {
        ids = sample(dim(X)[1], Nna_bycol, replace=FALSE)
        X[ids,j]=NA
    }
    return(X)
}

Adj_rand_index_table <- function(X, ncomponent, k,  miss_rate = seq(0.0, 0.4, 0.05), fo="./adj.csv", seed=1)
{
  ## Store results of different simulation
  Adjusted_Rand_index_table = data.frame( matrix(rep(0, 8*length(miss_rate)), ncol=8))

  row.names(Adjusted_Rand_index_table) = miss_rate
  names(Adjusted_Rand_index_table) = c("GMM_filtered", "GMM", "missForest", "KNN", "MICE","MICE_filtered", "mean", "median")

  for(mr in miss_rate){

      Mt2=NULL

      X_miss = X
      X_miss[,1:ncomponent] = insertNA(as.matrix(X_miss[, 1:ncomponent]), mr, seed)
      print(X_miss)
      print(head(X_miss))
      print(class(X_miss))
      Mt2 = compute_data_cluster_robust(as.data.frame(X_miss[,1:ncomponent]), k)
      X_miss["cluster"] = Mt2$clustered_data[,"cluster"]
      X_miss["entropy"] = Mt2$clustered_data[,"entropy"]
      if(!is.null(Mt2))
      {
        Adjusted_Rand_index_table[as.character(mr), "GMM"] = adjustedRandIndex(X_miss$component, X_miss$cluster)
        # Rand index when we filter out ambiguous Assignation
        well_identified = which(X_miss$entropy < 0.25)
        Adjusted_Rand_index_table[as.character(mr), "GMM_filtered"] = adjustedRandIndex(X_miss[well_identified, ]$component, X_miss[well_identified, ]$cluster)
      }

      #######################################################
      #######################################################
      # Evaluation basic imputation strategies (mean, median)
      #######################################################
      X_mean = X_miss[,1:ncomponent]
      X_median = X_miss[,1:ncomponent]

      X_mean[,1:ncomponent] = impute_by_fun(X_miss[,1:ncomponent])
      X_median[,1:ncomponent] = impute_by_fun(X_miss[,1:ncomponent], median)

      print(X_miss)

      Mt_mean <- compute_data_cluster_robust(X_mean[,1:ncomponent], k)
      Mt_median <- compute_data_cluster_robust(X_median[,1:ncomponent], k)
      if(!is.null(Mt_mean)){Adjusted_Rand_index_table[as.character(mr), 'mean'] = adjustedRandIndex(Mt_mean$clustered_data[,"cluster"], X$component)}
      if(!is.null(Mt_median)){Adjusted_Rand_index_table[as.character(mr), 'median'] = adjustedRandIndex(Mt_median$clustered_data[,"cluster"], X$component)}

      #######################################################
      #######################################################
      # Evaluation of MICE imputation
      # Let's trap the mice ;)
      #######################################################
      #######################################################

      X_imputed = X_miss[,1:ncomponent]
      nsamp = dim(X)[1]
      row.names(X_imputed) = 1:nsamp
      n_imp = 5
      imputation = mice(X_imputed[,1:ncomponent], m=n_imp, maxit = 50, method = 'pmm', seed = 500)
      # store cluster assignation for every imputation
      cl_mice = data.frame(matrix(0,dim(X)[1], ncol=n_imp))
      names(cl_mice) = 1:n_imp


      for(i in 1:n_imp){
          dat = complete(imputation, i)
          Mt3 = compute_data_cluster_robust(dat, k)
          if(!is.null(Mt3))
          {
            sort_clust = rank(sapply(Mt3$MNMmix@Means, function(x){x[[1]]}))
            cl_mice[,i] = sort_clust[Mt3$MNMmix@Assignments$Assignment]
          }else{
            print(paste0("Clustering failed for mice ",i))
          cl_mice[i] = NULL
          }
      }

      #X_imputed[,1:ncomponent] = complete(imputation)
      X_imputed["cluster"] = apply(cl_mice, 1, getmode)
      X_imputed["Proba_mode"] = apply(cl_mice[1:nsamp,]==X_imputed$cluster[1:nsamp], 1, mean)

      Adjusted_Rand_index_table[as.character(mr), 'MICE'] = adjustedRandIndex(X_imputed$cluster, X$component)
      Adjusted_Rand_index_table[as.character(mr), 'MICE_filtered'] = adjustedRandIndex(X_imputed[X_imputed$Proba_mode > 0.5,]$cluster, X[X_imputed$Proba_mode > 0.5,]$component)

      ################################################################
      # Evaluate missForest
      ################################################################
      if(ncomponent==2)
      {

        RF_imp <- missForest(cbind(X_miss[,1:ncomponent], rep(1,nsamp)))
      }
      else
      {
        RF_imp <- missForest(X_miss[,1:ncomponent])
      }
      X_RF = RF_imp$ximp[,1:ncomponent]
      Mt4 <- compute_data_cluster_robust(X_RF[,1:ncomponent], k)
      if(!is.null(Mt4)){
      X_RF["cluster"]= Mt4$MNMmix@Assignments$Assignment
      #X_RF["entropy"] = apply(Mt4@Responsibilities[,-1], 1, shannon_entropy)
      Adjusted_Rand_index_table[as.character(mr), 'missForest'] = adjustedRandIndex(X_RF$cluster, X$component)
    }

      ############################################################
      # kNN
      ############################################################
      X_kNN <- kNN(X_miss[,1:ncomponent])
      Mt5 <- compute_data_cluster_robust(X_kNN[, 1:ncomponent], k)

      X_kNN["cluster"]= Mt5$MNMmix@Assignments$Assignment
      #X_kNN["entropy"] = apply(Mt5@Responsibilities[,-1], 1, shannon_entropy)
      Adjusted_Rand_index_table[as.character(mr), 'KNN'] = adjustedRandIndex(X_kNN$cluster, X$component)

      print(Adjusted_Rand_index_table)

      write.csv(Adjusted_Rand_index_table, fo)
  }
    write.csv(Adjusted_Rand_index_table, fo)
}