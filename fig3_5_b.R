#The following code is the code used in Example 3.5.1 to and produces Figure 3.5 b)

#This code uses the package "ggplot2" written by Haldey Wickham and Winston Chang, published in 2016.
#The url for this package is https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf
#This code uses the package "GenSA" written by Sylvain Gubian, Yang Xiang, Brian Suomela and Julia Hoeng, published in 2018.
#The url for this package is https://cran.r-project.org/web/packages/GenSA/GenSA.pdf



library(GenSA)
library(ggplot2)

#Schweffel function
Schwef <- function(xx){
  d <- length(xx)
  sum <- sum(xx*sin(sqrt(abs(xx))))
  y<- 418.9829*d - sum
  return(y)
}
set.seed(1234)
dimension <- 30
global.min <- 0
tol <- 1e-5 #will stop searching with this absolute tolerance
lower <- rep(-500, dimension)
upper <- rep(500, dimension)

#run simulated annealing
sa1 <- GenSA(lower = lower, upper = upper, fn = Schwef, control=list(maxit=100))
sa1[c("value", "par", "counts")]

#Figure 3.5 b)
dat <- as.data.frame(sa1$trace.mat)
ggplot(dat, aes(x = nb.steps, y =current.minimum)) +
  geom_line(data = dat, 
            aes(x = nb.steps, y = function.value), 
            position = position_jitter(height = 0, width = 0.005),
            alpha = 1, colour ="red")+
  #geom_step(alpha = 0.3, colour="red") +
  theme(text = element_text(size=28)) +
  xlab("Iteration") +
  ylab("Minimum at Current Iteration")


