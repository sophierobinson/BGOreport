#The first section of the following code produces plots for Figure 1.1, Figure 2.1, Figure 2.2.
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

#squared exponential covariance function
calcSigma1 <- function(X1,X2,l=1) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
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

#toy objective function (Eqn. 1.2)
eg <- function(x){
  sin(x) + sin((2/3)*x)
}


#plotting the objective function (Figure 1.1)
xo <- seq(3.1, 20.4, len=50)
yo <- eg(xo)
plot(xo,yo)
x1 <- melt(xo)
y1 <- melt(yo)
obj <- data.frame(x = x1, y = y1)
objectivefunction <- ggplot(obj)+
  geom_line(aes(xo,yo),color="red")+
  theme_bw()
objectivefunction
#shows objective function and 3 functions drawn from GP using squared exponential function
a <- ggplot(valuesa,aes(x=x,y=value), obj) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  geom_line(aes(x1,y1), linetype=2) +
  theme(legend.position="none")+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 

#objective function as dotted line (for plots)
objfn <- ggplot(obj)+
  geom_line(aes(x1,y1), linetype=2)
objfn
#training set
x2 <- c(4, 10, 17)
y2 <- eg(x2)
f <- data.frame(x=x2, y=y2)
x <- f$x
k.xx <- calcSigma1(x,x)
k.xxs <- calcSigma1(x,x.star)
k.xsx <- calcSigma1(x.star,x)
k.xsxs <- calcSigma1(x.star,x.star)
sigma.n <- 0.00001

# Recalculate the mean and covariance functions
f.bar.star <- k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs

n.samples <- 5000
valuesb <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)

for (i in 1:n.samples) {
  valuesb[,i] <- mvrnorm(1, f.bar.star, cov.f.star)
}
valuesb <- cbind(x=x.star,as.data.frame(valuesb))
valuesb <- melt(valuesb,id="x")
#posterior distribution given training set with noise
b <- ggplot(valuesb ,aes(x=x,y=value), obj) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star,y=f.bar.star), aes(x=x,y=y),colour="#E69F00") +
  geom_point(data=f,aes(x=x,y=y)) +
  geom_line(aes(x1,y1), linetype=2) +
  theme_bw() +
  scale_y_continuous(lim=c(-4,4), name="output, f(x)") +
  xlab("input, x")

#Add a further point for evaluation x=13
x3 <- c(4, 10, 17, 13)
y3 <- eg(x3)
f <- data.frame(x=x3, y=y3)
x <- f$x
k.xx <- calcSigma1(x,x)
k.xxs <- calcSigma1(x,x.star)
k.xsx <- calcSigma1(x.star,x)
k.xsxs <- calcSigma1(x.star,x.star)
sigma.n <- 0.00001

# Recalculate the mean and covariance functions
f.bar.star2 <- k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
cov.f.star2 <- k.xsxs - k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs


n.samples <- 5000
valuesc <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)

for (i in 1:n.samples) {
  valuesc[,i] <- mvrnorm(1, f.bar.star2, cov.f.star2)
}
valuesc <- cbind(x=x.star,as.data.frame(valuesc))
valuesc <- melt(valuesc,id="x")
#posterior distribution given training set with further point and noise
c <- ggplot(valuesc ,aes(x=x,y=value), obj) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star,y=f.bar.star2), aes(x=x,y=y),colour="#E69F00") +
  geom_point(data=f,aes(x=x,y=y)) +
  geom_line(aes(x1,y1), linetype=2) +
  theme_bw() +
  scale_y_continuous(lim=c(-3,4), name="output, f(x)") +
  xlab("input, x")

#Figure 2.1
require(gridExtra)
grid.arrange(a, b, c)

#The following code is for the Acquisition functions (Chapter 3) and produces plots for Figure 3.1, Figure 3.2, Figure 3.3.

#GP-UCB ACQUSITION FUNCTION
UCB <- function(muNew, stdNew, t, d, v=1, delta=0.1){
  Kappa = sqrt(v*(2*log((t^(d/2 + 2))*(pi^2)/(3*delta))))
  return(muNew + Kappa*stdNew)
}
#variance of predictive distribution
std <- as.matrix(sqrt(diag(cov.f.star)))

#Figure 3.3
ucbacq0 <- UCB(f.bar.star, std, 10, 2, v=1, delta=0.1)
ucbacq1 <- melt(ucbacq0)
ucbacq2 <- data.frame(ucbacq1)
gg2 <- ggplot(ucbacq1) + 
  geom_line(aes(x=xo, y=ucbacq0), colour="blue")+
  geom_point(aes(x=12.985714, y=max(ucbacq0)), color="red")+
  theme_bw()
gg2
#finding maximum of the acquisiton
which.max(ucbacq0)
#29
#search xo for 29th element as x entry in geom_point
grid.arrange(b, gg2)


#finding fMax
currentbest <- max(y2)
#probability of improvement acqusition function
PI <- function(muNew, stdNew, fMax, epsilon){
  Z= (muNew - fMax - epsilon)/stdNew
  return(pnorm(Z))
}
#Figure 3.1
PIacq0 <- PI(f.bar.star, std, currentbest, 0.1)
PIacq1 <- melt(PIacq0)
PIacq2 <- data.frame(PIacq1)
gg3 <- ggplot(PIacq2)+
  geom_line(aes(x=xo, y=PIacq0), colour="blue")+
  geom_point(x=12.985714, y=max(PIacq0), colour="red")+
  theme_bw()
gg3
grid.arrange(b, gg3)
#finding maximum of PI
which.max(PIacq0)
#29

#Expected Improvement
EI <- function(muNew, stdNew, fMax, epsilon){
  Z= (muNew - fMax - epsilon)/stdNew
  fun = (muNew - fMax - epsilon)*pnorm(Z) + stdNew*dnorm(Z)
  return(fun)
}
#Figure 3.2
EIacq0 <- EI(f.bar.star, std, currentbest, 0.1)
EIacq1 <- melt(EIacq0)
EIacq2 <- data.frame(EIacq1)
gg4 <- ggplot(EIacq2)+
  geom_line(aes(x=xo, y=EIacq0), colour="blue")+
  geom_point(x=12.985714, y=max(EIacq0), colour="red")+
  theme_bw()
gg4
#finding max of EI for geom_point
which.max(EIacq0)
#29
grid.arrange(b, gg4)

