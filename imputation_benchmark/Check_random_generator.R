
setwd("/home/hjulienn/Project/missing_value_benchmark/imputation_benchmark/")
source("./Adjusted_rand_index_table.R")

D = data.frame(matrix(1,3,ncol=5))
dim(D)
for(i in 1:10){
  #print(insertNA(D,0.1,i))
  print(Adj_rand_index_table(D,5,6,seed=i))

}
