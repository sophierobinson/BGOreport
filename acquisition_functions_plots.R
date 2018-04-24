library(GPfit)
require(MASS)
require(gridExtra)
require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)
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

eg <- function(x){
  sin(x) + sin((2/3)*x)
}

x <- seq(3.1,20.4, len=50)
y <- eg(x)
plot(x,y)
x1 <- melt(x)
y1 <- melt(y)
obj <- data.frame(x = x1, y = y1)
objectivefunction <- ggplot(obj)+
  geom_line(aes(x,y),color="red")+
  theme(text = element_text(size=25))
objectivefunction

#shows objective function and 3 functions drawn from a GP prior distributions using squared exponential function
fig2_1 <- ggplot(valuesa,aes(x=x,y=value), obj) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable, color=variable)) +
  geom_line(aes(x1,y1), linetype=2) +
  theme(legend.position="none",text = element_text(size=25))+
  scale_y_continuous(lim=c(-3,3), name="f(x)") +
  xlab("input, x") 

fig2_1

#objective function as dotted line
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
sigma.n <- 0.1

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
fig2_1_b <- ggplot(valuesb ,aes(x=x,y=value), obj) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star,y=f.bar.star), aes(x=x,y=y),colour="#E69F00") +
  geom_point(data=f,aes(x=x,y=y)) +
  geom_line(aes(x1,y1), linetype=2) +
  theme(text = element_text(size=25)) +
  scale_y_continuous(lim=c(-4,4), name="f(x)") +
  xlab("input, x")
fig2_1_b


#GP-UCB ACQUSITION FUNCTION
x <- seq(3.1, 20.4, len=50)
UCB <- function(muNew, stdNew, t, d, v=1, delta=0.1){
  Kappa = sqrt(v*(2*log((t^(d/2 + 2))*(pi^2)/(3*delta))))
  return(muNew + Kappa*stdNew)
}
#variance of predictive distribution
std <- as.matrix(sqrt(diag(cov.f.star)))

ucbacq0 <- UCB(f.bar.star, std, 10, 2, v=1, delta=0.1)
ucbacq1 <- melt(ucbacq0)
ucbacq2 <- data.frame(ucbacq1)
fig3_3 <- ggplot(ucbacq1) + 
  geom_line(aes(x=x, y=ucbacq0), colour="blue")+
  geom_point(aes(x=12.985714, y=max(ucbacq0)), color="red")+
  theme(text = element_text(size=25))+
  ylab("GP-UCB")
fig3_3

#finding maximum of the acquisiton
which.max(ucbacq0)
#29
#search x for 29th element as x entry in geom_point
#Figure 3.3
require(gridExtra)
grid.arrange(fig2_1_b, fig3_4)



#finding fMax
currentbest <- max(y2)
#probability of improvement acqusition function
PI <- function(muNew, stdNew, fMax, epsilon){
  Z= (muNew - fMax - epsilon)/stdNew
  return(pnorm(Z))
}
PIacq0 <- PI(f.bar.star, std, currentbest, 0.2)
PIacq1 <- melt(PIacq0)
PIacq2 <- data.frame(PIacq1)
fig3_1_b <- ggplot(PIacq2)+
  geom_line(aes(x=x, y=PIacq0), colour="blue")+
  geom_point(x=12.985714, y=max(PIacq0), colour="red")+
  theme(text = element_text(size=25))+
  ylab("PI")

fig3_1_b

grid.arrange(fig2_1_b, fig3_1_b)

#finding maximum of PI
which.max(PIacq0)
#29

PIacq0 <- PI(f.bar.star, std, currentbest, 0)
PIacq1 <- melt(PIacq0)
PIacq2 <- data.frame(PIacq1)
fig3_1_a <- ggplot(PIacq2)+
  geom_line(aes(x=x, y=PIacq0), colour="blue")+
  geom_point(x=12.985714, y=max(PIacq0), colour="red")+
  theme(text = element_text(size=25))+
  ylab("PI")

fig3_1_a
grid.arrange(fig2_1_b, fig3_1_a)

#EXPECTED IMPROVEMENT
EI <- function(muNew, stdNew, fMax, epsilon){
  Z= (muNew - fMax - epsilon)/stdNew
  fun = (muNew - fMax - epsilon)*pnorm(Z) + stdNew*dnorm(Z)
  return(fun)
}
EIacq0 <- EI(f.bar.star, std, currentbest, 0.01)
EIacq1 <- melt(EIacq0)
EIacq2 <- data.frame(EIacq1)

fig3_2 <- ggplot(EIacq2)+
  geom_line(aes(x=x, y=EIacq0), colour="blue")+
  geom_point(x=12.985714, y=max(EIacq0), colour="red")+
  theme(text = element_text(size=25))+
  ylab("EI")

fig3_2

#finding max of EI for geom_point
which.max(EIacq0)
#29

grid.arrange(fig2_1_b, fig3_2)
