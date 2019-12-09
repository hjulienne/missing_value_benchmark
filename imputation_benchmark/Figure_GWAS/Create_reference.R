library(MGMM)
library(FactoMineR)
library("factoextra")
library(RNOmni)

shannon_entropy <- function(x){return(-sum(x*log2(x+10^-20)))}

setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_GWAS/")
source("../clustering_analysis_function.R")
Reg = read.csv("./raw_data/Regions_CARDIO_small_dim.csv")
Reg = Reg[1:1703,]
row.names(Reg) = Reg$snp_ids

Ndim = sum(grepl("z_", names(Reg)))

reg_w = which(Reg$UNIVARIATE_MIN_PVAL < 10^-8 | Reg$JASS_PVAL < 10^-8)
length(reg_w)
Reg = Reg[reg_w, grepl("z_", names(Reg))]
Reg = Reg[,1:Ndim]
Reg = align_zscores(Reg)

#Reg = as.data.frame(apply(Reg,2,rankNorm))
# Rand index when we filter out ambiguous Assignation
row.names(Reg) = Reg$snp_ids

S0=list()
Ndim = dim(Reg)[2]
for(i in 1:Nclust){
    S0[[i]] = matrix(rep(0,Ndim*Ndim), ncol=Ndim)
    diag(S0[[i]]) = 1
    # choose a random data point
}

Nclust=3
dim(Reg)
Cl = compute_data_cluster_robust(Reg, Nclust)

Reg["Cluster"] = Cl$clustered_data$cluster
Reg["entropy"] = Cl$MNMmix@Assignments[,"Entropy"]#Cl$clustered_data$entropy
entropy_treshold=0.1

Reg_pca = Reg[which(Reg$entropy < entropy_treshold),]

R = rank(tapply(Reg[,1], Reg$Cluster, mean))
Reg$Cluster = R[Reg$Cluster]

pca1 <- PCA(Reg_pca[,1:Ndim], scale.unit=TRUE, graph=F)

col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E", "#DD25B5", "#614ad3", "#23a13e", "#ed1330")

p1 = fviz_pca(pca1,
  geom=c("point"),
  pointsize=3,
  habillage= as.factor(Reg_pca$Cluster)) + scale_colour_manual(values=col_pal)

p2 = fviz_pca(pca1,
  geom=c("point"),
  pointsize=3, axes=c(2,3),
  habillage= as.factor(Reg_pca$Cluster)) + scale_colour_manual(values=col_pal)

plot_grid(p1,p2, nrow=1)

head(Reg)
ggplot(Reg[Reg$entropy<0.1,], aes(x=z_GLG_LDL, y=z_GIANT_BMI, color=as.factor(Cluster))) + geom_point()
ggplot(Reg[Reg$entropy<0.1,], aes(x=z_GLG_LDL, y=z_CARDIOGRAMPLUSC4D_CAD, color=as.factor(Cluster))) + geom_point()

write.csv(Reg[Reg$entropy<0.1,], "GWAS_cl.csv")
