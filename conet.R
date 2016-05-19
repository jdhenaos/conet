library(affy)
library(siggenes)
library(GEOquery)
library(vsn)

alzGSM <- read.table("Alzheimer_Chips.txt",stringsAsFactors = F)

############# FUNCION ##############

GetInfo <- function(GSE,GPL,dir="."){
  
  setwd(dir)
  
  for(i in GSE){
    getGEOSuppFiles(i)
  }
  
  files <- dir(".")[grep("GSE",dir("."))]
  
  for(j in files){
    untar(paste0(j,"/",j,"_RAW.tar"), exdir = paste0(j,"/"))
  }
  
  getGEOfile(GPL, destdir = ".")
}

####################################

############# FUNCION ##############

getaffy <- function(GSE){
  raw <- read.table(file = paste0(GSE,"/","filelist.txt"),sep = "\t",
                    header = T,comment.char = "#",stringsAsFactors = F)
  GSMs <- raw[,2]
  GSMs <- GSMs[grep(".CEL",GSMs)]
  affy <- ReadAffy(filenames = as.character(GSMs), compress = T,
                   celfile.path = GSE)
  
}


###################################

gene <- GeneSymbol("GPL570")

########## FUNCION ################

difexprs <- function(affy,treatment){
  rma <- rma(affy)
  eset <- ProbeFilter(rma,gene)
  y=treatment
  matrix <- as.matrix(eset)
  sam <- sam(matrix,y)
}


#################################

######### FUNCION ###############
plot(sam,1)
sum <- summary(sam,1,entrez=F)
dif.exp <- sum@row.sig.genes

################################

