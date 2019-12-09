library(MGMM)
library(colorspace)
library(RNOmni)

setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_GWAS/")
source("../clustering_analysis_function.R")

Reg = read.csv("./raw_data/Regions_CAD_risk.csv")
Reg = Reg[1:1703,]

names(Reg)

reg_w = which(Reg$UNIVARIATE_MIN_PVAL < 10^-6 | Reg$JASS_PVAL < 10^-6)

Reg = Reg[reg_w, grepl("z_", names(Reg))]

row.names(Reg) = Reg$snp_ids

Reg = align_zscores(Reg)
# Rand index when we filter

col_pal = c("#FF9D00", "#0C60BA","#008238","#EA401E", "#DD25B5")
col_heat = diverge_hsv(200)

get_trait <- function(x){strsplit(x, "_")[[1]][3]}
label_col <- names(Reg)[grepl("z_", names(Reg))]
label_col <- sapply(label_col, get_trait)
class(apply(Reg,2,rank))

heatmap.2(apply(Reg,2,rankNorm),
    margins=c(7, 4),
    trace="none",na.color="gray61",
    scale="none", key=TRUE,
    col = col_heat, Rowv=TRUE, Colv=TRUE, labRow=FALSE,labCol = label_col,
    cexCol=0.9, dendrogram="both",
    key.title = "Zscore", density.info="none",
   key.xlab = "Zscore",
   key.ylab = "")

D = dist(apply(Reg,2,rankNorm))
h1 = hclust(D)
plot(h1)
rect.hclust(h1 , k = 6, border = 2:6)
