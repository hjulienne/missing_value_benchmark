library(MGMM)
library(FactoMineR)
library("factoextra")
library("ggplot2")

Sigma = list()
means = list()
Nclust = 3
Sigma[[1]] <- matrix(c(2.5, 2.0, 0,
                        2.0, 2.5, 0,
                        0,0,0.3), 3, 3)

Sigma[[2]] <- matrix(c(4.4, 0.1, 0,
                        0.1, 0.3, 0,
                          0,0, 0.3),3,3)

Sigma[[3]] <- matrix(c(0.2, 0, 0,
                      0, 0.2, 0,
                      0,0,4.5),3,3)
norm(Sigma[[3]])

for(i in 1:Nclust){
  means[[i]] = c(0,0,0)
}

X = rGMM(900, d=3, k=3, M = means, S=Sigma)
component = row.names(X)
X = as.data.frame(X)
names(X) = c("y1", 'y2', 'y3')
X["component"] = component

cor_H0 = cor(X[,1:3])
inv_cor = solve(cor_H0)
chiZ = diag((as.matrix(X[,1:3]) %*% as.matrix(inv_cor)) %*%t(as.matrix(X[,1:3])))
p_val = 1-pchisq(chiZ,df=3)

p1 = ggplot(X[p_val<0.05,], aes(x=y1, y=y2, color=component)) + geom_point()
p2 = ggplot(X[p_val<0.05,], aes(x=y2, y=y3, color=component)) + geom_point()

write.csv(X[p_val<0.05,], "./Raw_data/data_GWAS_simulation.csv")
