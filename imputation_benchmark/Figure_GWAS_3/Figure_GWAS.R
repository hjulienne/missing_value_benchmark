library(colorspace)
library(gplots)
library(png)
setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_GWAS_3/")
source("../Imputation_perf_functions.R")

# Plot heatmap with reference clusters
dat2 = read.csv("GWAS_cl.csv", row.names=1)
dat2
#col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E", "#DD25B5", "#6417ff")
col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E", "#DD25B5", "#614ad3", "#23a13e", "#ed1330")


get_trait <- function(x){strsplit(x, "_")[[1]][3]}
label_col <- names(dat2)[grepl("z_", names(dat2))]
label_col <- sapply(label_col, get_trait)
lmat = matrix(c(3,0,5,1, 5, 2, 3,2, 4,2), ncol=5)


max_amp = max(abs(range(dat2[,grepl("z_",names(dat2))],na.rm=TRUE)))

# Create a linear color scale
max_tresh = 6
if(max_amp > max_tresh){
    bks = c(-max_amp, seq(-max_tresh, max_tresh, 0.5), max_amp)
}else{
    bks = c(-max_amp, seq(-max_tresh, max_tresh, 0.5), max_amp)
}

col_heat = diverge_hsv(length(bks)-1)


png("./heatmap_compressed_Z_score.png" , width=6, height=8, units = 'in', res = 300)

heatmap.2(as.matrix(dat2[order(dat2$Cluster,decreasing = TRUE), grepl("z_", names(dat2))]),
    #breaks = bks,
    margins=c(7, 4),
    trace="none",na.color="gray61",
    scale="none", key=TRUE,breaks=bks,
    col = col_heat, Rowv=FALSE, Colv=FALSE, labRow=FALSE,labCol = label_col,
    RowSideColors = col_pal[sort(dat2$Cluster,decreasing = TRUE)],
    cexCol=0.9, dendrogram="none",
    key.title = "Zscore", density.info="none",
   key.xlab = "Zscore",
   key.ylab = "", lmat=lmat, lwid=c(0.5,1,2.5,4,4), lhei=c(2,10))
dev.off()


png("./performance_GWAS.png" , width=6, height=8, units = 'in', res = 300)
#plot_performance_dashboard(data_perf)
pattern = "Adj_RI_GWAS_[0-9]+.csv"
plot_performance_adj_error("./benchmark_data/",pattern)
dev.off()

png("./fraction_data.png" , width=6, height=8, units = 'in', res = 300)
#plot_performance_dashboard(data_perf)
pattern = "Adj_RI_GWAS_[0-9]+.csv"

plot_performance_adj_error("./benchmark_data/",
pattern, methods=c("GMM_fraction_complete","MICE_fraction_complete","GMM_filtered_fraction","MICE_filtered_fraction"))

dev.off()

perf_name <- "./performance_GWAS.png"
heatmap_name <- "./heatmap_compressed_Z_score.png"
# example image
img1 <- readPNG(heatmap_name)
img2 <- readPNG(perf_name)
png("./Fig_panel_GWAS_2.png" , width=12, height=8, units = 'in', res = 300)
par(mar=rep(0,4)) # no margins
# layout the plots into a matrix w/ 12 columns, by row
layout(matrix(1:2, ncol=2, byrow=TRUE))

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(img1,0,0,1,1)

plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(img2,0,0,1,1)

dev.off()
