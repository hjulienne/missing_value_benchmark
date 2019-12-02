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

setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_simulation_GWAS/")
source("../Imputation_perf_functions.R")
source("../Adjusted_rand_index_table.R")

X = read.csv("./Raw_data/data_GWAS_simulation.csv")
X$X=NULL

ndim=3
ncomponent =3
names(X)[1:3] = c("z_1", "z_2", "z_3")
nsamp = dim(X)[1]


X[,1:ndim] = align_zscores(X[,1:ndim])
Xtest = insertNA(X, 0.1, 34)
# Reconstruction from complete data with classical GMM
col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E")

p1 = ggplot(X, aes(x=y1, y=y2, color=as.factor(component))) + geom_point() + geom_density2d() + theme_linedraw()+ lg_style
p1 = p1 + scale_colour_manual(values=col_pal)
p1 = p1 + lg_style + guides(color=guide_legend(title="gaussian component"))

X$component = as.factor(X$component)

p1 = ggplot(X, aes(x=z_1, y=z_2, color=component)) + geom_point() + scale_colour_manual(values=col_pal) + lg_style
p2 = ggplot(X, aes(x=z_2, y=z_3, color=component)) + geom_point()+ scale_colour_manual(values=col_pal) + lg_style

legend <- get_legend(p1 + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
p_row = plot_grid(p1 + theme(legend.position="none"),
                 p2+ theme(legend.position="none"))#+ lg_style + guides(color=guide_legend(title="gaussian component"))

png("Benchmark_simulation.png", width=10, height=5.5, unit="in", res=180)
  plot_grid(legend,p_row, nrow=2, align="v", rel_heights=c(0.1,0.9))
dev.off()

miss_rate = seq(0.0, 0.2, 0.025)
i=1

for(i in 21:40)
{
  fo = paste0("./benchmark_data/Adj_RI_simulation_GWAS_",i,".csv")
  SR = as.numeric(substr(as.character(as.numeric(Sys.time())), 13,17))
  Adj_rand_index_table(X,ncomponent,3 ,miss_rate, fo=fo, seed=SR)
}
