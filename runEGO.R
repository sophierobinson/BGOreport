#The following code was used to run Bayesian Global Optimization for groundwater contamination source localization (Chapter 4)
#It produces the Figures 4.3, 4.4, 4.5.
#This work is based on the code that can be found in the repository of GitHub at
#https://github.com/gpirot/BGICLP, by Pirot, G. and Krityakierne, T. and Ginsbourger, D. and Renard, P, published in 2017.

#This code uses the package "DiceKriging" published by Olivier Roustant, David Ginsbourger, Yves Deville. Contributors: Clement Chevalier,
#Yann Richet in 2015.
#The url for this package is https://cran.r-project.org/web/packages/DiceKriging/DiceKriging.pdf
#This code uses the package "DiceOptim" published by V. Picheny, D. Ginsbourger, O. Roustant
#with contributions by M. Binois, C. Chevalier, S. Marmin, and T. Wagner in 2016.
#The url for this package is https://cran.r-project.org/web/packages/DiceOptim/index.html
#This code uses the package "DiceDesign" published by Jessica Franco, Delphine Dupuy, Olivier Roustant, Guillaume
#Damblin and Bertrand Iooss in 2015.
#The url for this package is https://cran.r-project.org/web/packages/DiceDesign/DiceDesign.pdf
#This code uses the package "colorRamps" published by Tim Keitt in 2007.
#The url for this package is https://cran.r-project.org/web/packages/colorRamps/colorRamps.pdf

rm(list=ls())
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
#source("C:/Users/robo/Documents/Sophie/uni work/3rd year/project/BGICLP-master/src/functionGenerator.R")
#install.packages("methods")
#install.packages("DiceKriging")
#install.packages("DiceOptim")
#install.packages("DiceDesign")
#install.packages("lattice")
#install.packages("colorRamps")
library("methods")
library("DiceKriging")
library("DiceOptim")
library("lattice")
library("colorRamps")
source("C:/Users/robo/Documents/Sophie/uni work/3rd year/project/BGICLP-master/src/image.scale.R") # for plotting reasons only
set.seed(123)

# ***************************************************************
# Graphic parameters
# ***************************************************************
pal.1=colorRampPalette(c("yellow", "red", "black"), space="rgb")
colpoints <- "dodgerblue"#"cyan";#"blue";#"orange"; 
pchpoints <- 16
colnewpoint <- "orange"; 
colnewpointEI <- "red"; 
pchnewpoint <- 16
cexsize=1.5  


# ***************************************************************
# DEFINE SCENARIO AND EGO PARAMETERS
# ***************************************************************
# Geology: 1 or 2
selectedGeology <- 2
# Source c(89,-36) or c(100,10)
sourceCoord <- c(100,10)
# Norm: 1 for L1 norm, 2 for L2 norm
pNorm <- 1
# Wells ID: combination of wells identified individually between 1 and 25
selectedWells <- c(3,2,4,1,5,8,7,9,6,10,13,12,14,11,15,18,17,19,16,20,23,22,24,21,25)
# number of EGO iterations
nseq <- 100
# numberof points in the initial design
n0 <- 9
# seed for latin hypercube sampling design of initial points
lhsseed <- 1

# ***************************************************************
# BUILD OBJECTIVE FUNCTION 
# ***************************************************************

functionGenerator <- function(selectedWells,selectedGeology,sourceCoord,pNorm){
  geolName <- switch(selectedGeology,"A0","A4")
  file_name <- "C:/Users/robo/Documents/Sophie/uni work/3rd year/project/BGICLP-master/data/grid_25_wells_p1_A4_100_10.txt"
  
  out <- read.table(file_name, header=F)
  if (length(selectedWells)==1){
    out <- out[,selectedWells]^1/pNorm
  } else {
    out <- rowSums(out[,selectedWells])^1/pNorm
  }  
  return(out)
}


full_grid_par <- read.table("C:/Users/robo/Documents/Sophie/uni work/3rd year/project/BGICLP-master/data/full_grid_searching_zone_par.txt",header=F)
response.grid <- functionGenerator(selectedWells,selectedGeology,sourceCoord,pNorm)

argmin <- full_grid_par[which.min(response.grid),]
myTitle <- paste("Case 1")
#myTitle <- paste("Geology #",selectedGeology,", L",pNorm," norm, (x_s,y_s)=(",sourceCoord[1],",",sourceCoord[2],")",sep="")
xgrid <- seq(20,170,by=3)
ygrid <- seq(-75,75,by=3)
ngridx <- length(xgrid)
ngridy <- length(ygrid)

#pdf(file="EGOdemo1_final.pdf",width=12/2.54, height=11/2.54)
#layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(3,1))
library("DiceDesign")

	# GENERATE initial Latin Hypercube Sampling Design
 	design <- maximinESE_LHS(lhsDesign(n=n0, dimension=2, seed=lhsseed)$design)$design #round 
  design <- round(design*51+0.5)
  IDs <- (design[,2]-1)*51 + design[,1] 
   
   mydesign0 <- full_grid_par[IDs,]
   myresponse0 <- response.grid[IDs]
   mydesign <- mydesign0
   myresponse <- myresponse0
   minresponse <- rep(min(myresponse),nrow(mydesign0))
   mindistance <- rep(min(sqrt((mydesign0[,1]-argmin[,1])^2+ (mydesign0[,2]-argmin[,2])^2)),nrow(mydesign0))

   # RUN EGO and save figures in pdf (Figure 4.3)
   for(i in seq(1,nseq))
   {
      # GP model fitting and prediction 
      mykm <- km(~1, design=mydesign, response=myresponse,
              covtype="matern3_2", control=list(pop.size=50,trace=FALSE), 
              parinit=c(50, 50), lower=c(15,15), nugget = 1e-12)
      mypred <- predict(object=mykm, newdata=full_grid_par, type="UK")

      myEI <- apply(full_grid_par, 1, EI,mykm)
      EI.grid <- matrix(myEI, ngridx, ngridy)
      ################ SUMMARY GRAPHS ############
      kmean.grid <- matrix(mypred$mean, ngridx, ngridy) 
      breaks <- seq(min(matrix(kmean.grid, ngridx, ngridy)), max(matrix(kmean.grid, ngridx, ngridy)),length.out=100)
      par(mar=c(4.2,4.5,2.5,.5))
      image(xgrid, ygrid, matrix(kmean.grid, ngridx, ngridy),col=pal.1(length(breaks)-1),
            main=paste("GP mean prediction iter. ",i,sep=""),xlab="X(m)",ylab="Y(m)", axes=TRUE,cex.lab=cexsize,cex.axis=cexsize) # 
      points(mydesign[,1], mydesign[,2], pch=pchpoints, cex=1.5, col=colpoints)
      nextind <- which.max(myEI)
      box()
      par(mar=c(4.2,.5,2.5,5))
      image.scale(matrix(kmean.grid, ngridx, ngridy), col=pal.1(length(breaks)-1), breaks=breaks, horiz=FALSE,yaxt="n",xlab="",ylab="")
      axis(4,las=2,cex.axis=cexsize)
      box()
      contour(xgrid,ygrid,kmean.grid, nlevel=30, add=TRUE)
      #
      ksd.grid <- matrix(mypred$sd, ngridx, ngridy)
      breaks <- seq(min(matrix(ksd.grid, ngridx, ngridy)), max(matrix(ksd.grid, ngridx, ngridy)),length.out=100)
      par(mar=c(4.2,4.5,2.5,.5))
      image(xgrid, ygrid, matrix(ksd.grid, ngridx, ngridy),col=pal.1(length(breaks)-1),
      main=paste("GP standard deviation iter. ",i,sep=""),xlab="X(m)",ylab="Y(m)", axes=TRUE,cex.lab=cexsize,cex.axis=cexsize) # 
      box()
      par(mar=c(4.2,.5,2.5,5))
      image.scale(matrix(ksd.grid, ngridx, ngridy), col=pal.1(length(breaks)-1), breaks=breaks, horiz=FALSE,yaxt="n",xlab="",ylab="")
      axis(4,las=2,cex.axis=cexsize)
      box()
      contour(xgrid,ygrid,ksd.grid, nlevel=30, add=TRUE)
      points(mydesign[,1], mydesign[,2], pch=pchpoints, lwd=2, col=colpoints)
      #
      breaks <- seq(min(matrix(EI.grid, ngridx, ngridy)), max(matrix(EI.grid, ngridx, ngridy)),length.out=100)
      par(mar=c(4.2,4.5,2.5,.5))
      image(xgrid, ygrid, matrix(EI.grid, ngridx, ngridy),col=pal.1(length(breaks)-1),
            main=paste("Expected Improvement iter. ",i,sep=""), xlab="X(m)",ylab="Y(m)", axes=TRUE,cex.lab=cexsize,cex.axis=cexsize) #
      contour(xgrid,ygrid,EI.grid, nlevel=30, add=TRUE)
      points(mydesign[,1], mydesign[,2], pch=pchpoints, lwd=2, col=colpoints)
      nextind <- which.max(myEI)
      nextpoint <- full_grid_par[nextind,]
      points(nextpoint[,1], nextpoint[,2], pch=pchnewpoint, cex=1.5,col=colnewpointEI)
      box()
      par(mar=c(4.2,.5,2.5,5))
      image.scale(matrix(EI.grid, ngridx, ngridy), col=pal.1(length(breaks)-1), breaks=breaks, horiz=FALSE,yaxt="n",xlab="",ylab="")
      axis(4,las=2,cex.axis=cexsize)
      box()
      ###############################################################################
      # NEXT
      mydesign <- rbind(mydesign, nextpoint)
      myresponse <- c(myresponse, response.grid[nextind])
      minresponse <- c(minresponse,min(myresponse))
      mindistance <- c(mindistance,min(sqrt((mydesign[,1]-argmin[,1])^2+ (mydesign[,2]-argmin[,2])^2)))
   }
   dev.off()
   graphics.off()

   res.design <- mydesign[seq(1,nrow(mydesign0)+nseq),]
   res.cummin <- minresponse[seq(1,nrow(mydesign0)+nseq)]
   res.mindist <- mindistance[seq(1,nrow(mydesign0)+nseq)]
 #  SAVE DESIGN
   save(res.design,res.cummin,res.mindist,file=paste("resultsEGOlhs4_100",n0,"+",nseq,".Rda",sep=""))

#Figure 4.4
#for plotting final scenario with source and minimum overlayed, run 100 iterations and use the gpmean plot
hist(res.mindist)
kmean.grid <- matrix(mypred$mean, ngridx, ngridy)
breaks <- seq(min(matrix(kmean.grid, ngridx, ngridy)), max(matrix(kmean.grid, ngridx, ngridy)),length.out=100)
img1 <- {
par(mar=c(4.2,4.5,2.5,.5))
image(xgrid, ygrid, matrix(kmean.grid, ngridx, ngridy),col=pal.1(length(breaks)-1),
      main=paste("Case 4: Geology 2 and Source B"),xlab="X(m)",ylab="Y(m)", axes=TRUE,cex.lab=2,cex.axis=2, cex.main = 2) # 
points(mydesign[,1], mydesign[,2], pch=pchpoints, cex=1.5, col=colpoints)
points(100,10, pch =19, col = "black", cex=2, lwd=3)
points(74,15, pch = 4, col = "purple", cex=2, lwd=3)
nextind <- which.max(myEI)
box()
}
