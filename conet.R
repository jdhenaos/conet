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

difexprs <- function(affy,treatment,GPL){

  rma <- rma(affy)
  gene <- GeneSymbol(GPL)
  print("summarizing")
  eset <- ProbeFilter(rma,gene)
  matrix <- as.matrix(eset)
  pint("Differential analysis")
  sam <- sam(matrix,treatment)
  return(sam)
}


#################################

######### FUNCION ###############

getdifexprs <- function(sam,delta){
  plot(sam,delta)
  sum <- summary(sam,delta,entrez=F)
  dif.exp <- sum@row.sig.genes
  return(dif.exp)
}

################################

Aarray <- getaffy(GSE = "GSE8216")
t <- c(rep(1,6),rep(0,3))
Adife <- difexprs(affy = Parray,treatment = t)
Agenes <- getdifexprs(sam = Pdife,delta = 0.5)

###############################

GetInfo(GSE = "GSE8216",GPL = "GPL2025")
testarray <- getaffy(GSE = "GSE8216")
testtrait <- c(1,1,1,0,0,0)
testdif <- difexprs(affy = testarray,treatment = testtrait,GPL = "GPL2025")
