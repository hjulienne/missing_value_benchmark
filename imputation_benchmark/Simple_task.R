library(MASS)
library(ggplot2)
library(cowplot)
library(mice)
library(mclust)
library(MNMix)
library(MGMM)
library("missForest")
library('VIM')

#source('../clustering_analysis_function.R')

setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/")
source("./Imputation_perf_functions.R")

Nsamp = 2000

###
# Create a bivariate dataset from a gaussian distribution
# This is just a simple task to benchmark different imputation methods
###

shannon_entropy <- function(x){return(-sum(x*log2(x)))}
mu1 = c(0,0)


means = list(mu1)
Sigma <- matrix(c(5, 4, 4, 5),2,2)

X = rGMM(Nsamp, d=2, k=1, M = means, S=Sigma)

X = as.data.frame(X)
names(X) = c("y1", 'y2')


# Reconstruction from complete data with classical GMM
save_plot( "simple_task.pdf", ggplot(X, aes(x=y1, y=y2)) + geom_point(color="royalblue") + geom_density2d(color="orangered") + lg_style
)
insertNA <- function(X, fraction)
{
    Nna_bycol = ceiling(fraction * dim(X)[1])

    for(j in 1:dim(X)[2])
    {
        ids = sample(dim(X)[1], Nna_bycol, replace=FALSE)
        X[ids,j]=NA
    }
    return(X)
}

miss_rate = seq(0.05, 0.5, 0.05)
## Store results of different simulation
dist_table = data.frame( matrix(rep(0, 4*length(miss_rate)), ncol=4))

row.names(dist_table) = miss_rate
names(dist_table) = c("missForest", "KNN", "MICE","MICE_filtered")



for(mr in miss_rate){

  X_miss = X[,1:2]
  X_miss[,1:2] = insertNA(as.matrix(X_miss[,1:2]), mr)


    #######################################################
    #######################################################
    # Evaluation of MICE imputation
    #######################################################
    #######################################################

    X_imputed = X_miss[,1:2]
    dim(X_imputed)
    row.names(X_imputed) = 1:Nsamp
    names(X_imputed)

    n_imp = 5
    imputation = mice(X_imputed[,1:2], m=n_imp, maxit = 50, method = 'pmm', seed = 500)

    # store cluster assignation for every imputation
    N_missing = sum(is.na(X_imputed))
    cl_mice = data.frame(matrix(0, N_missing, ncol=n_imp))
    cor_p = rep(0, n_imp)

    for(i in 1:n_imp){
        dat = complete(imputation, i)
        cor_p[i] = cor(dat[is.na(X_miss)], X[is.na(X_miss)])
        cl_mice[,i] = dat[is.na(X_miss)]
    }

    X_imputed[,1:2] = complete(imputation)

    im_sd = apply(cl_mice, 1, sd)
    Qc = quantile(im_sd, p=0.5)
    imp_val = apply(cl_mice, 1, mean)
    true_val = X[is.na(X_miss)]
    D_imp = data.frame(imputed_value = imp_val, real_value=true_val)

    dist_table[as.character(mr), 'MICE'] = cor(imp_val, true_val)
    dist_table[as.character(mr), 'MICE_filtered'] = cor(imp_val[im_sd < Qc], true_val[im_sd < Qc])


    save_plot(paste0("./simple_task/MICE_", mr, ".pdf"),
    ggplot(D_imp, aes(x=imputed_value, y=real_value)) + geom_point() + geom_density2d()
  )

    save_plot(paste0("./simple_task/MICE_filtered_", mr, ".pdf"),
    ggplot(D_imp[im_sd < Qc,], aes(x=imputed_value, y=real_value)) + geom_point() + geom_density2d()
  )

    ################################################################
    # Evaluate missForest
    ################################################################

    RF_imp <- missForest(cbind(X_miss[,1:2], rep(1, dim(X_miss)[1])))
    X_RF = RF_imp$ximp[,1:2]

    imp_val = X_RF[is.na(X_miss)]
    true_val = X[is.na(X_miss)]
    D_imp = data.frame(imputed_value = imp_val, real_value=true_val)

    dist_table[as.character(mr), 'missForest'] = cor(imp_val, true_val, use='pairwise.complete.obs')

    save_plot(paste0("./simple_task/missForest_", mr, ".pdf"),
    ggplot(D_imp, aes(x=imputed_value, y=real_value)) + geom_point() + geom_density2d()
  )

    ############################################################
    # kNN
    ############################################################
    X_kNN <- kNN(X_miss[,1:2])


    imp_val = X_kNN[,1:2][is.na(X_miss)]
    true_val = X[is.na(X_miss)]
    D_imp = data.frame(imputed_value = imp_val, real_value=true_val)
    length(imp_val)
    dist_table[as.character(mr), 'KNN'] = cor(imp_val, true_val)

    save_plot(paste0("./simple_task/KNN_", mr, ".pdf"),
    ggplot(D_imp, aes(x=imputed_value, y=real_value)) + geom_point() + geom_density2d()
  )

}


write.csv(dist_table, 'correlation_simple_task.csv')
