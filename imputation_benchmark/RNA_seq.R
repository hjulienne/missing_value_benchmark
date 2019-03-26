library(MASS)
library(ggplot2)
library(cowplot)
library(mice)
library(mclust)
library(MNMix)
library(MGMM)
library("missForest")
library('VIM')
library(FactoMineR)
library("factoextra")
#source('../clustering_analysis_function.R')
setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/")

source("./Imputation_perf_functions.R")

###
# Create a bivariate dataset from a gaussian mixture with 4 components
# The component are placed so knowledge of the first variables should be enough
# to retrieve while the second variable alone shouldn't be enough
# The Covariance matrix is chosen so the local correlation inside one component
# goes in a different direction than the global correlation
###
shannon_entropy <- function(x){return(-sum(x*log2(x)))}

X = read.csv("../datasets/Other_examples/TCGA-PANCAN/data_top_20.csv", row.names=1)
tumour = read.csv("../datasets/Other_examples/TCGA-PANCAN/labels.csv", row.names=1)

nsamp = dim(X)[1]
nvar = sum(grepl("gene_", names(X)))
#sort(apply(X[grepl("gene_", names(X))], 2, sd), decreasing)

X["component"] = tumour[,1]
unique(tumour[,1])

pca1 <- PCA(X[, grepl("gene_", names(X))], scale.unit=TRUE, graph=T)

Mt = fit.MNMix(as.matrix(X[, grepl("gene_", names(X))]), k=5, parallel=FALSE, report=TRUE, maxit=500)
Mt_bis = fit.GMM(as.matrix(X[, grepl("gene_", names(X))]), k=5, parallel=FALSE, report=TRUE, maxit=500)

X["cluster"]= Mt@Assignments$Assignment
X["cluster_new"]= Mt_bis@Assignments$Assignment
X["entropy"] = apply(Mt@Responsibilities[,-1], 1, shannon_entropy)

p1 = fviz_pca_ind(pca1,
  geom=c("point"),
  pointsize=3,
  habillage= as.factor(tumour[,1]))

p2 = fviz_pca_ind(pca1,
  geom=c("point"),axes=c(2,3),
  pointsize=3,
  habillage= as.factor(tumour[,1]))

# Reconstruction from complete data with classical GMM

pdf("rna_seq.pdf", width=12, height=5)
plot_grid(p1,p2)
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

miss_rate = seq(0.05, 0.3, 0.05)
## Store results of different simulation
Adjusted_Rand_index_table = data.frame( matrix(rep(0, 6*length(miss_rate)), ncol=6))

row.names(Adjusted_Rand_index_table) = miss_rate
names(Adjusted_Rand_index_table) = c("GMM_filtered", "GMM", "missForest", "KNN", "MICE","MICE_filtered")

entropy_treshold = 0.25
mr = 0.05



for(mr in miss_rate){
  tryCatch(
    {
    X_miss = X
    ids_miss = sample(dim(X)[1], dim(X)[1]*0.25)
    X_miss[ids_miss,grepl("gene_", names(X))] = insertNA(as.matrix(X_miss[ids_miss, grepl("gene_", names(X))]), ((mr * 4.0)/ 3.0))

    Mt2 = fit.MNMix(as.matrix(X_miss[,grepl("gene_", names(X))]),k=4)
    X_miss["cluster"]= Mt2@Assignments$Assignment
    X_miss["entropy"] = apply(Mt2@Responsibilities[,-1], 1, shannon_entropy)

    Adjusted_Rand_index_table[as.character(mr), "GMM"] = adjustedRandIndex(X_miss$component, X_miss$cluster)

    # Rand index when we filter out ambiguous Assignation
    Adjusted_Rand_index_table[as.character(mr), "GMM_filtered"] = adjustedRandIndex(X_miss[X_miss$entropy < 0.25, ]$component, X_miss[X_miss$entropy < 0.25, ]$cluster)

    #######################################################
    #######################################################
    # Evaluation of MICE imputation
    # Let's trap the mice ;)
    #######################################################
    #######################################################


    X_imputed = X_miss[,1:nvar]
    row.names(X_imputed) = 1:nsamp

    n_imp = 5
    imputation = mice(X_imputed[,grepl("gene_", names(X))], m=n_imp, maxit = 50, method = 'pmm', seed = 500)


    # store cluster assignation for every imputation
    cl_mice = data.frame(matrix(0,dim(X)[1], ncol=n_imp))

    for(i in 1:n_imp){
        dat = complete(imputation, i)
        Mt3 = fit.MNMix(as.matrix(dat), k=4)

        sort_clust = rank(sapply(Mt3@Means, function(x){x[[1]]}))

        cl_mice[,i] = sort_clust[Mt3@Assignments$Assignment]
    }

    X_imputed[,grepl("gene_", names(X))] = complete(imputation)
    head(X_imputed)
    X_imputed["cluster"] = apply(cl_mice, 1, getmode)
    dim(cl_mice)
    dim(X_imputed)

    X_imputed["Proba_mode"] = apply(cl_mice[1:nsamp,]==X_imputed$cluster[1:nsamp], 1, mean)

    Adjusted_Rand_index_table[as.character(mr), 'MICE'] = adjustedRandIndex(X_imputed$cluster, X$component)
    Adjusted_Rand_index_table[as.character(mr), 'MICE_filtered'] = adjustedRandIndex(X_imputed[X_imputed$Proba_mode > 0.5,]$cluster, X[X_imputed$Proba_mode > 0.5,]$component)

    ################################################################
    # Evaluate missForest
    ################################################################

    RF_imp <- missForest(cbind(X_miss[,grepl("gene_", names(X))], rep(1, dim(X_miss)[1])))
    X_RF = RF_imp$ximp[,grepl("gene_", names(X))]

    head(X_RF)

    Mt4 <- fit.MNMix(as.matrix(X_RF[grepl("gene_", names(X))]), k=4)

    X_RF["cluster"]= Mt4@Assignments$Assignment
    X_RF["entropy"] = apply(Mt4@Responsibilities[,-1], 1, shannon_entropy)
    Adjusted_Rand_index_table[as.character(mr), 'missForest'] = adjustedRandIndex(X_RF$cluster, X$component)



    ############################################################
    # kNN
    ############################################################


    X_kNN <- kNN(X_miss[,grepl("gene_", names(X))])

    Mt5 <- fit.MNMix(as.matrix(X_kNN[, 1:nvar]), k=4)

    X_kNN["cluster"]= Mt5@Assignments$Assignment
    X_kNN["entropy"] = apply(Mt5@Responsibilities[,-1], 1, shannon_entropy)
    Adjusted_Rand_index_table[as.character(mr), 'KNN'] = adjustedRandIndex(X_kNN$cluster, X$component)

    print(Adjusted_Rand_index_table)
  },   error=function(cond) {
            message(paste("Error for missing rate:", mr))

            message(cond)
            # Choose a return value in case of error
            return(NA)})
}


write.csv(Adjusted_Rand_index_table, 'Adj_rand_index_RNA_seq.csv')
