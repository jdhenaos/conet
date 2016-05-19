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
  
  files <- dir(".")[grep("^GSE[0-9]",dir("."))]
  
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
  return(affy)
}

###################################

########## FUNCION ################

difexprs <- function(affy,treatment){
  
  #gene <- GeneSymbol(GPL)
  rma <- rma(affy)
  print("summarizing")
  eset <- ProbeFilter(rma,gene)
  matrix <- as.matrix(eset)
  print("Differential analysis")
  sam <- sam(matrix,treatment)
  
  plot(sam,delta)
  sum <- summary(sam,delta,entrez=F)
  dif.exp <- sum@row.sig.genes
  return(sam)
}


#################################


################################

Aarray <- getaffy(GSE = "GSE28379")
t <- c(rep(0,8),rep(1,22))
#gene <- GeneSymbol("GPL570")
Adife <- difexprs(affy = Aarray,treatment = t)
Agenes <- getdifexprs(sam = Adife,delta = 1)

###############################





testarray <- getaffy(GSE = "GSE8216")
testtrait <- c(1,1,1,1,0,0,0,0)
testdif <- difexprs(affy = testarray,treatment = testtrait,GPL = "GPL570")
testgenes <- getdifexprs(sam = testdif,delta = 1)
