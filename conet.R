library(affy)
library(siggenes)
library(GEOquery)
library(vsn)

alzGSM <- read.table("Descargas/Alzheimer_Chips.txt",stringsAsFactors = F)

for(i in alzGSM$V1){
  getGEOSuppFiles(i)
}