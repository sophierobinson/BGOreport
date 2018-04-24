library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#matern1/2 with noise var = 1
calcSigma11 <- function(X1, X2, l=1){
  Sigma11 <- matrix(rep(0, length(X1)*length(X2)), nrow=(length(X1)))
  for(i in 1:nrow(Sigma11)){
    for(j in 1:ncol(Sigma11)){
      r = abs(X1[i]-X2[j])
      Sigma11[i,j] <- exp(-r/l)
    }
  }
  return(Sigma11)
}

x.star11 <- seq(-5,5, len=50)
sigma11 <- calcSigma11(x.star11, x.star11)
n.samples <- 3
valuesa1 <- matrix(rep(0,length(x.star11)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa1[,i] <- mvrnorm(1, rep(0, length(x.star11)), sigma11)
}
valuesa1 <- cbind(x=x.star11,as.data.frame(valuesa1))
valuesa1 <- melt(valuesa1,id="x")

#Figure 2.3 a)
maternprior1 <- ggplot(valuesa1,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-1.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none",text = element_text(size=30))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 
maternprior1

