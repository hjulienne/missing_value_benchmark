library(MGMM)
library(RNOmni)

setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_GWAS/")
source("../clustering_analysis_function.R")

Reg = read.csv("./raw_data/Regions_CARDIO_small_dim.csv")
Reg = Reg[1:1703,]

reg_w = which(Reg$UNIVARIATE_MIN_PVAL < 10^-8 | Reg$JASS_PVAL < 10^-8)

Reg = Reg[reg_w, grepl("z_", names(Reg))]
row.names(Reg) = Reg$snp_ids

Reg = align_zscores(Reg)
Reg

#Reg = as.data.frame(apply(Reg,2,rankNorm))
# Rand index when we filter out ambiguous Assignation
plot(Reg$z_GIANT_BMI, Reg$z_GIANT_WHR)
perf_data = compute_perf_by_k_bootstrap(Reg,
     maxclust=15,
    repetition=10)

write.csv(perf_data,"performance_data.csv")

perf_data = read.csv("./performance_data.csv",stringsAsFactors=FALSE)

tapply(perf_data$Silhouette, perf_data$k, mean)
png(paste0("./clust_perf.png" ), width=1200, height=400)
    print(optimal_k_plot_bootstrap(as.data.table(perf_data[,-1])))
dev.off()
