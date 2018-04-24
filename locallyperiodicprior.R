library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#locally periodic covariance function
calcSigma6 <- function(X1, X2, l=2, p=0.5){
  Sigma6 <- matrix(rep(0, length(X1)*length(X2)), nrow=(length(X1)))
  for(i in 1:nrow(Sigma6)){
    for(j in 1:ncol(Sigma6)){
      r = abs(X1[i]-X2[j])
      a = sin(pi*abs(X1[i]-X2[j])/p)^2
      b <- exp((-2*a)/l^2)
      c <- exp(-(r^2)/(2*l^2))
      Sigma6[i,j] <- b*c
    }
  }
  return(Sigma6)
}


x.star6 <- seq(-10,20, len=50)
sigma6 <- calcSigma6(x.star6, x.star6)
n.samples <- 3
valuesa6 <- matrix(rep(0,length(x.star6)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa6[,i] <- mvrnorm(1, rep(0, length(x.star6)), sigma6)
}
valuesa6 <- cbind(x=x.star6,as.data.frame(valuesa6))
valuesa6 <- melt(valuesa6,id="x")

#Figure 2.5
localperiodicprior <- ggplot(valuesa6,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-1.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none",text = element_text(size=30))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 

localperiodicprior

