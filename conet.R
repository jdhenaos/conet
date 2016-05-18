library(affy)
library(siggenes)
library(GEOquery)
library(vsn)

alzGSM <- read.table("Descargas/Alzheimer_Chips.txt",stringsAsFactors = F)

############# FUNCION ##############

setwd("alzheimer/")

for(i in alzGSM$V1){
  getGEOSuppFiles(i)
}

files <- dir(".")[grep("GSE",dir("."))]

for(j in files){
  untar(paste0(j,"/",j,"_RAW.tar"), exdir = paste0(j,"/"))
}

getGEOfile("GPL570", destdir = ".")

####################################

############# FUNCION ##############

