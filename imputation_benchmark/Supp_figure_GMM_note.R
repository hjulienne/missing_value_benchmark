library(RColorBrewer)
library(data.table)
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

X <- read.csv("./data_bivariate_simulation.csv")#,row.names=1)
X = X[,-1]

# Reconstruction from complete data with classical GMM
col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E")

head(X)
p1 = ggplot(X, aes(x=y1, y=y2, color=as.factor(component))) + geom_point() + geom_density2d() + theme_linedraw()+ lg_style
p1 = p1 + scale_colour_manual(values=col_pal)
p1 = p1 + lg_style + guides(color=guide_legend(title="gaussian component"))

dat = read.csv("./Adj_rand_index.csv")
dat
dat_long= melt(dat, id.vars= "X")

p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line(lwd=2) + geom_point(size=2.5)
p = p + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_linedraw()
p = p + lg_style+ scale_colour_brewer(palette = "Set1")# brewer.pal(6,"Set1")
p = p + guides(color=guide_legend(title="method"))

png("Benchmark_simulation.png", width=10, height=5.5, unit="in", res=180)
  plot_grid(p1, p)
dev.off()
