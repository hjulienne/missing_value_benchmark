library(MASS)
library(ggplot2)
library(mice)

source('../clustering_analysis_function.R')
###

mu1 = c(1,1)
mu2 = c(4,1)

mu3 = c(-4,-1)
mu4 = c(-1, -1)
means = list(mu1, mu2, mu3, mu4)

Sigma <- matrix(c(0.9,-0.6,-0.6,0.9),2,2)

Nsamp = 400
X = rGMM(1000, d=2, k=4, M = means, S=Sigma )
X = as.data.frame(X)
X["component"] = row.names(X)
ggplot(as.data.frame(X), aes(x=y1, y=y2, color=component)) + geom_point() + geom_density2d()

#ggplot(data_cl, aes(x=X1, y=X2,color=cluster)) + geom_point() + geom_density2d()

Mt = fit.MNMix(as.matrix(data_cl[, grepl("X", names(data_cl))]), k=4, parallel=FALSE, report=TRUE, maxit=500)

data_cl["cluster"]= Mt@Assignments$Assignment
data_cl["entropy"] = apply(Mt@Responsibilities[,-1], 1, shannon_entropy)

# Perfect recomstruction
ggplot(data_cl, aes(x=X1, y=X2,color=class)) + geom_point() + geom_density2d()
ggplot(data_cl, aes(x=X1, y=X2,color=as.factor(cluster))) + geom_point() + geom_density2d()


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

# Let's trap the mice:
data_cl_miss = data_cl[,1:3]
data_cl_miss[,1:2] = insertNA(as.matrix(data_cl_miss[,1:2]), 0.4)

Mt2 = fit.MNMix(as.matrix(data_cl_miss[,1:2]),k=4)

data_cl_miss["cluster"]= Mt2@Assignments$Assignment
data_cl_miss["entropy"] = apply(Mt2@Responsibilities[,-1], 1, shannon_entropy)

table(data_cl_miss$class, data_cl_miss$cluster)

# Confusion matrix when we filter out ambiguous Assignation
table(data_cl_miss[data_cl_miss$entropy < 0.4, ]$class, data_cl_miss[data_cl_miss$entropy < 0.4, ]$cluster)


ggplot(data_cl, aes(x=X1, y=X2,color=class)) + geom_point() + geom_density2d()
ggplot(data_cl, aes(x=X1, y=X2,color=as.factor(cluster))) + geom_point() + geom_density2d()


data_cl_imputed = data_cl_miss[,1:3]

imputation = mice(data_cl_miss[,1:2], m=5, maxit = 50, method = 'pmm', seed = 500)
Mt3 = fit.MNMix(as.matrix(complete(imputation,1)),k=4)

data_cl_imputed["cluster"]= Mt3@Assignments$Assignment
data_cl_imputed["entropy"] = apply(Mt3@Responsibilities[,-1], 1, shannon_entropy)

table(data_cl_imputed$class, data_cl_imputed$cluster)
