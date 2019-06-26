library(data.table)
library(FactoMineR)
library("factoextra")
library("tsne")
library(ggplot2)
library(infotheo)
library(gplots)
library(mclust)
library(heatmap3)

data <- fread("~/Project/missing_value_benchmark/datasets/Other_examples/TCGA-PANCAN/data.csv")
label <- fread("~/Project/missing_value_benchmark/datasets/Other_examples/TCGA-PANCAN/labels.csv", header=TRUE)

dim(data)
table(label$Class)

# Let's plot random genes
n_gene = 10
id_rd = sample((dim(data)[2] -1), n_gene)

data_long = melt(data[, id_rd, with=FALSE])
ggplot(data_long, aes(x=value)) + geom_histogram() + facet_wrap( ~ variable)

gene_completion = apply(data[, 2:20532, with=FALSE][,.SD>0], 2, mean)

data.frame(gene_completion)
ggplot(data.frame(gene_completion), aes(x=gene_completion)) + geom_histogram() + geom_vline(xintercept = 0.05)

filled_gene = which(gene_completion > 0.1)
pca1 <- PCA(data[, filled_gene, with=FALSE], scale.unit=TRUE, graph=T)
pca1$eig[,3] < 75
fviz_pca_ind(pca1,
  geom=c("point"),
  pointsize=3,
  habillage= as.factor(label$Class))

fviz_pca_ind(pca1,
  geom=c("point"), axes=c(2,3),
  pointsize=3,
  habillage= as.factor(label$Class))

gene_pval <- rep(0, length(filled_gene))
i=1

get_pval <- function(x)
{
  pval <- anova(lm(x ~ label$Class))[1,5]
}
gene_pval = data[, sapply(.SD,get_pval), .SDcols = names(data)[filled_gene]]
# heatmap of the best 20 markers:
heatmap3(cor(data[,.SD, .SDcols = names(filled_gene[order(gene_pval)[1:20]])]),
  balanceColor=T)

pca2 <- PCA(data[,.SD, .SDcols = names(filled_gene[order(gene_pval)[1:20]])], scale.unit=TRUE, graph=T)

fviz_screeplot(pca2)
fviz_pca_ind(pca2,
    geom=c("point"),
    pointsize=3,
    habillage= as.factor(label$Class))

fviz_pca_ind(pca2,
    geom=c("point"), axes=c(2,3),
    pointsize=3,
    habillage= as.factor(label$Class))


write.csv(data[, .SD, .SDcols = names(filled_gene[order(gene_pval)[1:20]])], "~/Project/missing_value_benchmark/datasets/Other_examples/TCGA-PANCAN/data_top_20.csv" )
