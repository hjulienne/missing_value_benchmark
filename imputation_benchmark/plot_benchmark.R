library(ggplot2)
library(reshape)

dat = read.csv("./Adj_rand_index.csv")
dat_long= melt(dat, id.vars= "X")


 p = ggplot(dat_long, aes(x=X, y = value, color=variable)) + geom_line() + geom_point() + xlab("miss ratio") + ylab("Adjusted Rand Index") + theme_bw()
ggsave("benchmark.pdf", p)
