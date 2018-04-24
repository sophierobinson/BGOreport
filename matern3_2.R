library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#matern3/2 with noise var = 1
calcSigma2 <- function(X1,X2,l=1) {
  Sigma2 <- matrix(rep(0, length(X1)*length(X2)), nrow=(length(X1)))
  for(i in 1:nrow(Sigma2)){
    for(j in 1:ncol(Sigma2)){
      Sigma2[i,j] <- (1+(sqrt(3)*abs(X1[i]-X2[j])/l))*exp(-(sqrt(3)*abs(X1[i]-X2[j]))/l)
    }
  }
  return(Sigma2)
}

x.star2 <- seq(-5,5, len=50)
sigma2 <- calcSigma2(x.star2, x.star2)
n.samples <- 3
valuesa2 <- matrix(rep(0,length(x.star2)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa2[,i] <- mvrnorm(1, rep(0, length(x.star2)), sigma2)
}
valuesa2 <- cbind(x=x.star2,as.data.frame(valuesa2))
valuesa2 <- melt(valuesa2,id="x")

#Figure 2.3 b)
maternprior2 <- ggplot(valuesa2,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-1.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none")+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 
maternprior2

