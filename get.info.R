# get.info
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

get.info <- function(GSE,GPL,dir="."){
  
  # Moves to the specifed directory
  
  setwd(dir)
  
  # Download the raw data from GEO DataSets database
  
  sapply(as.vector(t(GSE)), getGEOSuppFiles)
  
  # Obtain the name of the compress data
  
  files <- dir(".")[grep("^GSE[0-9]",dir("."),ignore.case = T)]
  
  # Decompress all the raw data
  
  for(j in files){
    untar(paste0(j,"/",j,"_RAW.tar"), exdir = paste0(j,"/"))
  }
  
  # Download the .soft file of the microaaray chip
  
  getGEOfile(GPL, destdir = ".")
}