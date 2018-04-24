library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#periodic covariance function
calcSigma5 <- function(X1, X2, l=2, p=0.5){
  Sigma5 <- matrix(rep(0, length(X1)*length(X2)), nrow=(length(X1)))
  for(i in 1:nrow(Sigma5)){
    for(j in 1:ncol(Sigma5)){
      a = sin(pi*abs(X1[i]-X2[j])/p)^2
      b <- exp((-2*a)/l^2)
      Sigma5[i,j] <- b
    }
  }
  return(Sigma5)
}


x.star5 <- seq(-8,8, len=50)
sigma5 <- calcSigma5(x.star5, x.star5)
n.samples <- 3
valuesa5 <- matrix(rep(0,length(x.star5)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa5[,i] <- mvrnorm(1, rep(0, length(x.star5)), sigma5)
}
valuesa5 <- cbind(x=x.star5,as.data.frame(valuesa5))
valuesa5 <- melt(valuesa5,id="x")

#Figure 2.4
periodicprior <- ggplot(valuesa5,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-1.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none",text = element_text(size=35))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 
periodicprior
