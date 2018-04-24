#for case 1 (geology 1 and source A) the source was found in 53 iterations
#for case 2 (geology 1 and source B) the next closest point to the source is (101,9) this is found by iertion 34
#for case 3 (geology 2 and source A) the source was found n 18 iterations
#for case 4 (geology 2 and source B) the minimum found by the algorithm was in iteration 63 (74,15) this is 26.5 m away from the contaminant source


hist(res.mindist[1:100],
     main="After 100 iterations", 
     xlab="Distance to the source (m)", 
     ylab ="Counts",
     col="purple",
     xlim = c(0,35),
     cex.axis = 2.5,
     cex.lab = 2.5,
     cex.main = 2.5)
