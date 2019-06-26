library(MASS)
library(ggplot2)
library(cowplot)
library(mice)
library(mclust)
library(MGMM)
library("missForest")
library('VIM')
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

shannon_entropy <- function(x){return(-sum(x*log2(x+10^-20)))}

X <- read.csv("./data_bivariate_simulation.csv")#,row.names=1)
X = X[,-1]
Mt = fit.GMM(as.matrix(X[, grepl("y", names(X))]), k=4, parallel=FALSE, report=TRUE, maxit=100)

X["cluster"]= Mt@Assignments$Assignment

Mt@Responsibilities
apply(Mt@Responsibilities, 1, shannon_entropy)[1:5]

X["entropy"] = apply(Mt@Responsibilities, 1, shannon_entropy)
head(X)
# Reconstruction from complete data with classical GMM
col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E")

head(X)
p1 = ggplot(X, aes(x=y1, y=y2, color=as.factor(component))) + geom_point() + geom_density2d() + theme_linedraw()+ lg_style
p1 = p1 + scale_colour_manual(values=col_pal)
p1

p2 = ggplot(X, aes(x=y1, y=y2, color=as.factor(cluster))) + geom_point() + geom_density2d() + theme_linedraw()+ lg_style
p2 = p2 + scale_colour_manual(values=col_pal)
p3 = ggplot(X, aes(x=y1, y=y2, color= entropy)) + geom_point() + geom_density2d() + theme_linedraw()+ lg_style

pdf("clustering_gaussian_mixture.pdf", width=12, height=4)
plot_grid(p1,p2, p3, nrow=1)
dev.off()
print( paste0("GMM on the complete data RA ", adjustedRandIndex(X$component, X$cluster)))

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

miss_rate = seq(0.05, 0.4, 0.05)

## Store results of different simulation
Adjusted_Rand_index_table = data.frame( matrix(rep(0, 8*length(miss_rate)), ncol=8))

row.names(Adjusted_Rand_index_table) = miss_rate
names(Adjusted_Rand_index_table) = c("GMM_filtered", "GMM", "missForest", "KNN", "MICE","MICE_filtered", "mean", "median")

entropy_treshold = 0.25
mr =0.15

for(mr in miss_rate){

    X_miss = X[,1:3]
    X_miss[,1:2] = insertNA(as.matrix(X_miss[,1:2]), mr)

    Mt2 = fit.GMM(as.matrix(X_miss[,1:2]),k=4, maxit=20)
    X_miss["cluster"]= Mt2@Assignments$Assignment
    X_miss["entropy"] = apply(Mt2@Responsibilities[,-1], 1, shannon_entropy)

    Adjusted_Rand_index_table[as.character(mr), "GMM"] = adjustedRandIndex(X_miss$component, X_miss$cluster)


    save_plot(paste0("./scatter_plots/GMM_", mr, ".pdf"),
    plot_against_reference(X, X_miss)
  )
    # Rand index when we filter out ambiguous Assignation
    Adjusted_Rand_index_table[as.character(mr), "GMM_filtered"] = adjustedRandIndex(X_miss[X_miss$entropy < 0.25, ]$component, X_miss[X_miss$entropy < 0.25, ]$cluster)

    save_plot(paste0("./scatter_plots/GMM_filtered", mr, ".pdf"),
    plot_against_reference(X[X_miss$entropy < 0.25, ], X_miss[X_miss$entropy < 0.25, ])
  )
    #######################################################
    #######################################################
    # Evaluation basic imputation (mean, median)
    #
    #######################################################
    X_mean = X_miss[,1:3]
    X_median = X_miss[,1:3]
    head(X_miss)
  #  median(X_miss[,1], na.rm=TRUE)

    X_mean[,1:2] = impute_by_fun(X_miss[,1:2])
    X_median[,1:2] = impute_by_fun(X_miss[,1:2], median)

    Mt_mean <- fit.GMM(as.matrix(X_mean[1:2]), k=4,maxit=20)
    Mt_median <- fit.GMM(as.matrix(X_median[1:2]), k=4,maxit=20)

    Adjusted_Rand_index_table[as.character(mr), 'mean'] = adjustedRandIndex(Mt_mean@Assignments$Assignment, X$component)
    Adjusted_Rand_index_table[as.character(mr), 'median'] = adjustedRandIndex(Mt_median@Assignments$Assignment, X$component)

    #######################################################
    #######################################################
    # Evaluation of MICE imputation
    # Let's trap the mice ;)
    #######################################################
    #######################################################

    X_imputed = X_miss[,1:3]
    row.names(X_imputed) = 1:2000

    n_imp = 5
    imputation = mice(X_imputed[,1:2], m=n_imp, maxit = 50, method = 'pmm', seed = 500)

    # store cluster assignation for every imputation
    cl_mice = data.frame(matrix(0,dim(X)[1], ncol=n_imp))
    for(i in 1:n_imp){
        dat = complete(imputation, i)
        Mt3 = fit.GMM(as.matrix(dat), k=4,maxit=20)
        print(Mt3@Means)
        sort_clust = rank(sapply(Mt3@Means, function(x){x[[1]]}))
        print(sort_clust)
        cl_mice[,i] = sort_clust[Mt3@Assignments$Assignment]
    }

    X_imputed[,1:2] = complete(imputation)

    X_imputed["cluster"] = apply(cl_mice, 1, getmode)
    X_imputed["Proba_mode"] = apply(cl_mice[1:n_imp,]==X_imputed$cluster[1:n_imp], 1, mean)

    Adjusted_Rand_index_table[as.character(mr), 'MICE'] = adjustedRandIndex(X_imputed$cluster, X_imputed$component)
    Adjusted_Rand_index_table[as.character(mr), 'MICE_filtered'] = adjustedRandIndex(X_imputed[X_imputed$Proba_mode > 0.5,]$cluster, X_imputed[X_imputed$Proba_mode > 0.5,]$component)


    save_plot(paste0("./scatter_plots/MICE_", mr, ".pdf"),
    plot_reconstruction_against_reference(X, X_imputed)
  )

    save_plot(paste0("./scatter_plots/MICE_filtered_", mr, ".pdf"),
    plot_reconstruction_against_reference(X[X_imputed$Proba_mode > 0.5, ], X_imputed[X_imputed$Proba_mode > 0.5, ])
  )

    ################################################################
    # Evaluate missForest
    ################################################################

    RF_imp <- missForest(cbind(X_miss[,1:2], rep(1, dim(X_miss)[1])))
    X_RF = RF_imp$ximp[,1:2]

    Mt4 <- fit.GMM(as.matrix(X_RF[1:2]), k=4,maxit=20)

    X_RF["cluster"]= Mt4@Assignments$Assignment
    X_RF["entropy"] = apply(Mt4@Responsibilities[,-1], 1, shannon_entropy)
    Adjusted_Rand_index_table[as.character(mr), 'missForest'] = adjustedRandIndex(X_RF$cluster, X$component)

    save_plot(paste0("./scatter_plots/RF_", mr, ".pdf"),
    plot_reconstruction_against_reference(X, X_RF)
  )

    ############################################################
    # kNN
    ############################################################


    X_kNN <- kNN(X_miss[,1:2])
    Mt5 <- fit.GMM(as.matrix(X_kNN[1:2]), k=4,maxit=20)
    X_kNN["cluster"]= Mt5@Assignments$Assignment
    X_kNN["entropy"] = apply(Mt5@Responsibilities[,-1], 1, shannon_entropy)
    Adjusted_Rand_index_table[as.character(mr), 'KNN'] = adjustedRandIndex(X_kNN$cluster, X$component)

    save_plot(paste0("./scatter_plots/kNN_", mr, ".pdf"),
    plot_reconstruction_against_reference(X, X_RF)
  )

}

Adjusted_Rand_index_table

write.csv(Adjusted_Rand_index_table, 'Adj_rand_index.csv')
