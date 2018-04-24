library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#Vary l and sig2 in thie squared exponential function to produce different prior distributions
calcSigma1 <- function(X1,X2,l=3, sig2=1) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- sig2*exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}

x.star <- seq(3.1,20.4, len=50)
sigma <- calcSigma1(x.star, x.star)
n.samples <- 3
valuesa <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa[,i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
valuesa <- cbind(x=x.star,as.data.frame(valuesa))
valuesa <- melt(valuesa,id="x")

#Figure 2.6
fig2_6 <- ggplot(valuesa,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none", text = element_text(size=40))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 
fig2_6

