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

raw <- read.table(file = paste0("GSE16759/","filelist.txt"),sep = "\t",
                  header = T,comment.char = "#",stringsAsFactors = F)
GSMs <- raw[,2]
GSMs <- GSMs[grep(".CEL",GSMs)]
affy <- ReadAffy(filenames = as.character(GSMs), compress = T,celfile.path = "GSE16759")

