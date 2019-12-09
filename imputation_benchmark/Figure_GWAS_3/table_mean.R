library(data.table)



setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/Figure_GWAS_3/")

DT = fread("./performance_data.csv")


DT[,.(count=.N), by=k]
