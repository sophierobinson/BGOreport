functionGenerator <- function(selectedWells,selectedGeology,sourceCoord,pNorm){
  geolName <- switch(selectedGeology,"A0","A4")
  file_name <- "C:/Users/robo/Documents/Sophie/uni work/3rd year/project/BGICLP-master/data/grid_25_wells_p1_A0_89_-36.txt"
  
  out <- read.table(file_name, header=F)
  if (length(selectedWells)==1){
    out <- out[,selectedWells]^1/pNorm
  } else {
	out <- rowSums(out[,selectedWells])^1/pNorm
  }  
  return(out)
}
