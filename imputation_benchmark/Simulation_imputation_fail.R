library(MASS)
library(ggplot2)
library(cowplot)
library(mice)
library(DMwR)
library(mclust)
library(MNMix)
library("missForest")
library('VIM')

library(DMwR)
source('../clustering_analysis_function.R')
source("./Imputation_perf_functions.R")

###
# Create a bivariate dataset from a gaussian mixture with 4 components
# The component are placed so knowledge of the first variables should be enough
# to retrieve while the second variable alone shouldn't be enough
# The Covariance matrix is chosen so the local correlation inside one component
# goes in a different direction than the global correlation
###
mu1 = c(-4,-1)
mu2 = c(-1, -1)

mu3 = c(1,1)
mu4 = c(4,1)

means = list(mu1, mu2, mu3, mu4)
Sigma <- matrix(c(0.9, -0.6, -0.6, 0.9),2,2)

X = rGMM(2000, d=2, k=4, M = means, S=Sigma)
X = as.data.frame(X)
names(X) = c("y1", 'y2')
X["component"] = row.names(X)


Mt = fit.MNMix(as.matrix(X[, grepl("y", names(X))]), k=4, parallel=FALSE, report=TRUE, maxit=500)
X["cluster"]= Mt@Assignments$Assignment
X["entropy"] = apply(Mt@Responsibilities[,-1], 1, shannon_entropy)



# Reconstruction from complete data with classical GMM
p1 = ggplot(X, aes(x=y1, y=y2,color=component)) + geom_point() + geom_density2d() + lg_style
p2 = ggplot(X, aes(x=y1, y=y2,color=as.factor(cluster))) + geom_point() + geom_density2d() + lg_style
p3 = ggplot(X, aes(x=y1, y=y2,color= entropy)) + geom_point() + geom_density2d() + lg_style
plot_grid(p1,p2, p3)

print( paste0("GMM on the complete data RA", adjustedRandIndex(X$component, X$cluster)))

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

hist(X$entropy)

miss_rate = seq(0.05, 0.3, 0.05)
## Store results of different simulation
Adjusted_Rand_index_table = data.frame( matrix(rep(0, 6*length(miss_rate)), ncol=6))

row.names(Adjusted_Rand_index_table) = miss_rate
names(Adjusted_Rand_index_table) = c("GMM_filtered", "GMM", "missForest", "KNN", "MICE","MICE_filtered")

entropy_treshold = 0.25


for(mr in miss_rate){
    X_miss = X[,1:3]
    X_miss[,1:2] = insertNA(as.matrix(X_miss[,1:2]), mr)

    Mt2 = fit.MNMix(as.matrix(X_miss[,1:2]),k=4)
    X_miss["cluster"]= Mt2@Assignments$Assignment
    X_miss["entropy"] = apply(Mt2@Responsibilities[,-1], 1, shannon_entropy)

    Adjusted_Rand_index_table[as.character(mr), "GMM"] = adjustedRandIndex(X_miss$component, X_miss$cluster)

    pdf(paste0("./scatter_plots/GMM_", mr, ".pdf"))
    plot_against_reference(X, X_miss)
    dev.off()

    # Rand index when we filter out ambiguous Assignation
    Adjusted_Rand_index_table[as.character(mr), "GMM_filtered"] = adjustedRandIndex(X_miss[X_miss$entropy < 0.25, ]$component, X_miss[X_miss$entropy < 0.25, ]$cluster)

    pdf(paste0("./scatter_plots/GMM_filtered", mr, ".pdf"))
    plot_against_reference(X[X_miss$entropy < 0.25, ], X_miss[X_miss$entropy < 0.25, ])
    dev.off()

    #######################################################
    #######################################################
    # Evaluation of MICE imputation
    # Let's trap the mice ;)
    #######################################################
    #######################################################

    X_imputed = X_miss[,1:3]
    row.names(X_imputed) = 1:2000
    names(X_imputed)

    n_imp = 5
    imputation = mice(X_imputed[,1:2], m=n_imp, maxit = 50, method = 'pmm', seed = 500)

    # store cluster assignation for every imputation
    cl_mice = data.frame(matrix(0,dim(X)[1], ncol=n_imp))
    for(i in 1:n_imp){
        dat = complete(imputation, i)
        Mt3 = fit.MNMix(as.matrix(dat), k=4)
        print(Mt3@Means)
        sort_clust = rank(sapply(Mt3@Means, function(x){x[[1]]}))
        print(sort_clust)
        cl_mice[,i] = sort_clust[Mt3@Assignments$Assignment]
    }



    X_imputed["cluster"] = apply(cl_mice, 1, getmode)
    X_imputed["Proba_mode"] = apply(cl_mice[1:n_imp,]==X_imputed$cluster[1:n_imp], 1, mean)

    Adjusted_Rand_index_table[as.character(mr), 'MICE'] = adjustedRandIndex(X_imputed$cluster, X_imputed$component)
    Adjusted_Rand_index_table[as.character(mr), 'MICE_filtered'] = adjustedRandIndex(X_imputed[X_imputed$Proba_mode > 0.5,]$cluster, X_imputed[X_imputed$Proba_mode > 0.5,]$component)


    pdf(paste0("./scatter_plots/MICE_", mr, ".pdf"))
    plot_against_reference(X, X_imputed)
    dev.off()

    pdf(paste0("./scatter_plots/MICE_filtered_", mr, ".pdf"))
    plot_against_reference(X[X_imputed$Proba_mode > 0.5, ], X_imputed[X_imputed$Proba_mode > 0.5, ])
    dev.off()

    ################################################################
    # Evaluate missForest
    ################################################################

    RF_imp <- missForest(cbind(X_miss[,1:2], rep(1, dim(X_miss)[1])))
    X_RF = RF_imp$ximp[,1:2]

    Mt4 <- fit.MNMix(as.matrix(X_RF[1:2]), k=4)

    X_RF["cluster"]= Mt4@Assignments$Assignment
    X_RF["entropy"] = apply(Mt4@Responsibilities[,-1], 1, shannon_entropy)
    Adjusted_Rand_index_table[as.character(mr), 'missForest'] = adjustedRandIndex(X_RF$cluster, X$component)

    pdf(paste0("./scatter_plots/RF_", mr, ".pdf"))
    plot_against_reference(X, X_RF)
    dev.off()

    ############################################################
    # kNN
    ############################################################


    X_kNN <- kNN(X_miss[,1:2])
    Mt5 <- fit.MNMix(as.matrix(X_kNN[1:2]), k=4)
    X_kNN["cluster"]= Mt5@Assignments$Assignment
    X_kNN["entropy"] = apply(Mt5@Responsibilities[,-1], 1, shannon_entropy)
    Adjusted_Rand_index_table[as.character(mr), 'KNN'] = adjustedRandIndex(X_kNN$cluster, X$component)

    pdf(paste0("./scatter_plots/kNN_", mr, ".pdf"))
    plot_against_reference(X, X_RF)
    dev.off()
}
