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


