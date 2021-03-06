library(colorspace)
library(gplots)
library(png)
setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_simulation_GWAS/")
source("../Imputation_perf_functions.R")
source("../Adjusted_rand_index_table.R")
# Plot heatmap with reference clusters
dat2 = read.csv("./Raw_data/data_GWAS_simulation.csv")

dim(dat2)[1]
dat2$X = NULL
#col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E", "#DD25B5", "#6417ff")
col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E", "#DD25B5", "#614ad3", "#23a13e", "#ed1330")

col_heat = diverge_hsv(200)
names(dat2)[1:3] = c("z_1", "z_2", "z_3")
get_trait <- function(x){strsplit(x, "_")[[1]][3]}
label_col <- names(dat2)[grepl("z_", names(dat2))]
label_col <- sapply(label_col, get_trait)
lmat = matrix(c(3,0,5,1, 5, 2, 3,2, 4,2), ncol=5)
dat2[,1:3] = align_zscores(dat2[,1:3])

png("./heatmap_compressed_Z_score.png" , width=6, height=8, units = 'in', res = 500)
head(dat2)
heatmap.2(as.matrix(dat2[order(dat2$component,decreasing = TRUE), grepl("z_", names(dat2))]),
    #breaks = bks,
    margins=c(7, 4),
    trace="none",na.color="gray61",
    scale="none", key=TRUE,
    col = col_heat, Rowv=FALSE, Colv=FALSE, labRow=TRUE,
    RowSideColors = col_pal[sort(dat2$component,decreasing = TRUE)],
    cexCol=0.9, dendrogram="none",
    key.title = "Zscore", density.info="none",
   key.xlab = "Zscore",
   key.ylab = "", lmat=lmat, lwid=c(0.5,1,2.5,4,4), lhei=c(2,10))
dev.off()


png("./performance_GWAS_simu.png" , width=6, height=8, units = 'in', res = 500)
#plot_performance_dashboard(data_perf)
pattern = "Adj_RI_simulation_GWAS_[0-9]+.csv"

plot_performance_adj_error("./benchmark_data/",pattern)

dev.off()

png("./fraction_data.png" , width=6, height=8, units = 'in', res = 500)
#plot_performance_dashboard(data_perf)
pattern = "Adj_RI_simulation_GWAS_[0-9]+.csv"

plot_performance_adj_error("./benchmark_data/",
pattern, methods=c("GMM_fraction_complete","MICE_fraction_complete","GMM_filtered_fraction","MICE_filtered_fraction"))

dev.off()

perf_name <- "./performance_GWAS_simu.png"
heatmap_name <- "./heatmap_compressed_Z_score.png"
# example image
img1 <- readPNG(heatmap_name)
img2 <- readPNG(perf_name)
png("./panel_GWAS_simulation.png" , width=12, height=8, units = 'in', res = 500)
par(mar=rep(0,4)) # no margins
# layout the plots into a matrix w/ 12 columns, by row
layout(matrix(1:2, ncol=2, byrow=TRUE))

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(img1,0,0,1,1)

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(img2,0,0,1,1)

dev.off()
