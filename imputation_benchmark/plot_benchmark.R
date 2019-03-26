library(ggplot2)
library(reshape)


dat = read.csv("./Adj_rand_index.csv")
dat_long= melt(dat, id.vars= "X")

p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line() + geom_point() + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_bw()
ggsave("benchmark.png", p,width=5, height=4)

dat = read.csv("./Adj_rand_index_RNA_seq.csv")
dat_long= melt(dat, id.vars= "X")

p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line() + geom_point() + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_bw()
ggsave("benchmark_rna_seq.png", p,width=5, height=4)

dat = read.csv("./correlation_simple_task.csv")

dat_long= melt(dat, id.vars= "X")

 p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line() + geom_point() + xlab("missing data ratio") + ylab("Adjusted Rand Index") + theme_bw()

ggsave("benchmark_simple_task.png", p, width=5, height=4)


dat = read.csv("./correlation_moderate_task.csv")
dat
dat_long= melt(dat, id.vars= "X")
p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line() + geom_point() + xlab("miss ratio") + ylab("Adjusted Rand Index") + theme_bw()

ggsave("benchmark_moderate_task.png", p, width=5, height=4)
