#The following code was used to produce Figure 2.3 d)


#The GP regression section follows the code produced by James Keirstead: "Gaussian process regression with R" April 5, 2012.
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
