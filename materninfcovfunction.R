library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)

#materninf with noise var = 1
calcSigma4 <- function(X1, X2, l=1){
  Sigma4 <- matrix(rep(0, length(X1)*length(X2)), nrow=(length(X1)))
  for(i in 1:nrow(Sigma4)){
    for(j in 1:ncol(Sigma4)){
      r = abs(X1[i]-X2[j])
      Sigma4[i,j] <- exp(-(r^2)/(2*l^2))
    }
  }
  return(Sigma4)
}

x.star4 <- seq(-5,5, len=50)
sigma4 <- calcSigma4(x.star4, x.star4)
n.samples <- 3
valuesa4 <- matrix(rep(0,length(x.star4)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  valuesa4[,i] <- mvrnorm(1, rep(0, length(x.star4)), sigma4)
}
valuesa4 <- cbind(x=x.star4,as.data.frame(valuesa4))
valuesa4 <- melt(valuesa4,id="x")

#Figure 2.3 d)
maternprior4 <- ggplot(valuesa4,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-1.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  theme(legend.position="none",text = element_text(size=30))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 
maternprior4
