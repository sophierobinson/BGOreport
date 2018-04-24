
target <- function(X){ 
  t = 0.3*exp(-0.2*X^2) 
  return(t)
}
#plots the target distribution
x00 <- seq(-5,10,len=50)
y00 <- target(x00)
plot(x00,y00)
x <- melt(x00)
y <- melt(y00)
targ <- data.frame(x = x, y = y)

#Figure 2.7 a)
unimodal <- ggplot(targ)+
  geom_line(aes(x,y),color="green")+
  theme(text = element_text(size=40))
unimodal 
#proposal distribution
proposal <- function(Y,X){
  return(dnorm(X, Y, 100))
}
#generates Y
genprop <- function(X){
  return(rnorm(1,X,100))
}

#metropolis-hastings
mh = function(Nsim){
  X = rep(runif(1), Nsim)
  for(i in 2:Nsim){
    currentstate = X[i-1]
    Y = genprop(currentstate)
    acceptprob = (target(Y)*proposal(Y,currentstate))/(target(currentstate)*proposal(currentstate,Y))
    if(runif(1) < acceptprob){
      X[i] = Y
    } else {
      X[i] = currentstate
    }
  }
 return(X)
}

eg <- mh(100000)
par(mfrow=c(1,1))

#Figure 2.7 b)
plot(eg, type="l", xlab="Iteration", ylab="X", col=3, main=NULL, cex.lab=1.5, cex.axis=2)   #plots the iterations against approximations
#Figure 2.7 c)
hist(eg, xlab="MCMC samples", ylab = "frequency", col=3, main = NULL, cex.lab=1.6, cex.axis=2)  #histogram of approximations

