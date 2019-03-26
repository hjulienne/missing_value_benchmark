library(ggplot2)
library(MGMM)

Nsamp = 10

###
# Create a bivariate dataset from a gaussian distribution
# This is just a simple task to benchmark different imputation methods
###

shannon_entropy <- function(x){return(-sum(x*log2(x)))}
mu1 = c(0,0)


means = list(mu1)
Sigma <- matrix(c(5, 4, 4, 5),2,2)

cov2cor(Sigma)
X = rGMM(Nsamp, d=2, k=1, M = means, S=Sigma)
X = as.data.frame(X)
names(X) = c("y1", 'y2')



lm1 <- lm(X$y2 ~ X$y1, data = X)
X["predicted"] = predict(lm1)

X_missings = as.data.frame(rGMM(1, d=2, k=1, M = means, S=Sigma))
names(X_missings) = c("y1")

X_missings['predict_new'] = lm1$coefficients[1] + lm1$coefficients[2] * X_missings[,1]



p1 <- ggplot(X, aes(x = y1, y = y2)) +
  geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
  geom_segment(aes(xend = y1, yend = predicted), alpha = .2) +  # alpha to fade lines
  geom_point(size=2) +
  geom_point(aes(y = predicted), shape = 1) +
  theme_bw()  # Add theme for cleaner look
X_missings$y1

p1 = p1 + geom_point(data=X_missings, aes(x = y1), y=0, color="deepskyblue1", size=2) #+ geom_segment(aes(xend = y1, yend = predicted), alpha = .2)
p1 = p1 + geom_point(data=X_missings, aes(x = y1, y=predict_new), color="orangered1", size=2)
p1 = p1 + geom_segment(data=X_missings, aes(x=y1, y=0, xend = y1, yend = predict_new), alpha = .2)

p1

ggsave('PMM_illustration.png', p1,width = 4, height = 4)
