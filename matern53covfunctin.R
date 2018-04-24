library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#matern5/2 with noise var = 1
calcSigma3 <- function(X1, X2, l=1){
  Sigma3 <- matrix(rep(0, length(X1)*length(X2)), nrow=(length(X1)))
  for(i in 1:nrow(Sigma3)){
    for(j in 1:ncol(Sigma3)){
      r = abs(X1[i]-X2[j])
      Sigma3[i,j] <- (1+(sqrt(5)*r/l)+(5*r^2/3*l^2))*exp(-(sqrt(5)*r)/l)
    }
  }
  return(Sigma3)
}

x.star3 <- seq(-5,5, len=50)
sigma3 <- calcSigma3(x.star3, x.star3)
n.samples <- 3
valuesa3 <- matrix(rep(0,length(x.star3)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa3[,i] <- mvrnorm(1, rep(0, length(x.star3)), sigma3)
}
valuesa3 <- cbind(x=x.star3,as.data.frame(valuesa3))
valuesa3 <- melt(valuesa3,id="x")

#Figure 2.3 c)
maternprior3 <- ggplot(valuesa3,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-1.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none",text = element_text(size=30))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 
maternprior3
