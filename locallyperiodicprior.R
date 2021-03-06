#This code produces Figure 2.5, the Locally Periodic covariance function.

#The GP regression section of this code followa the code produced by James Keirstead: "Gaussian process regression with R" April 5, 2012.
#The code can be sourced here https://www.r-bloggers.com/gaussian-process-regression-with-r/

#This code uses the package "GPfit" written by Blake MacDonald and Pritam Ranjan and Hugh Chipman, published in 2015.
#The url is https://cran.r-project.org/web/packages/GPfit/index.html
#This code uses the package "MASS" written by Brian Ripley et. al., published in 2016.
#The url is https://cran.r-project.org/web/packages/MASS/MASS.pdf
#This code uses the package "gridExtra" written by Baptiste Augie, published in 2017.
#The url is https://cran.r-project.org/web/packages/gridExtra/index.html
#This code uses the package "plyr" written by Hadley Wickham, published in 2016.
#The url for this package is https://cran.r-project.org/web/packages/plyr/plyr.pdf
#This code uses the package "reshape2" written by Hadley Wickham, published in 2017.
#The url for this package is https://cran.r-project.org/web/packages/reshape2/reshape2.pdf
#This code uses the package "ggplot2" written by Haldey Wickham and Winston Chang, published in 2016.
#The url for this package is https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf


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

