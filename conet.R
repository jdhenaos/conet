library(affy)
library(siggenes)
library(GEOquery)
library(vsn)
library(igraph)

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

difexprs <- function(affy,treatment,fdr){
  #gene <- GeneSymbol(GPL)
  rma <- rma(affy)
  print("summarizing")
  #eset <- ProbeFilter(rma,gene)
  eset <- meanProbe(gene,rma)
  matrix <- as.matrix(eset)
  print("Differential analysis")
  sam <- sam(matrix,treatment)
  tab <- show(sam)
  mtab <- as.matrix(data.frame(tab$Delta,tab$FDR))
  filt <- mtab[mtab[,2]<=fdr,]
  delta <- as.numeric(filt[1,1])
  plot(sam,delta)
  sum <- summary(sam,delta,entrez=F)
  dife <- sum@row.sig.genes
  genes <- eset[dife,]
  return(genes)
}

#############################

####### FUNCION #############

CreateNet <- function(difexp){
  
  simil <- abs(cor(t(difexp),use =  "pairwise.complete.obs"))
  
  pcv <- seq(0.01,0.99,by = 0.01)
  
  Cis <- vector()
  C0s <- vector()
  
  count <- 1
  
  for (val in pcv) {
    ady <- matrix(0,ncol = ncol(simil), nrow = nrow(simil))
    
    for(i in 1:nrow(simil)){
      ady[which(simil[,i]>=val),i]<-1
      ady[which(simil[,i]<val),i]<-0
    }
    
    G = graph.adjacency(ady,mode="undirected",diag=FALSE)
    
    Ci <- transitivity(G,type = "globalundirected")
    
    if(is.nan(Ci)){Ci <- 0}
    
    K1 <- sum(degree(G,loops = F))
    K2 <- sum(degree(G,loops = F)^2)
    k1 <- (1/length(V(G)))*K1
    k2 <- (1/length(V(G)))*K2
    
    C0 <- ((k2-k1)^2)/(k1^3*length(V(G)))
    
    if(is.nan(C0)){C0 <- 0}
    
    Cis[count] <- Ci
    C0s[count] <- C0
    count <- count + 1
  }
  
  C <- vector()
  
  for (counter in 1:(length(pcv)-1)) {
    if((Cis[counter] - C0s[counter]) > (Cis[counter+1] - C0s[counter+1])){
      C[counter] <- pcv[counter]
    } 
  }
  
  C <- na.omit(C)
  
  fit <-1
  
  for (z in as.vector(C)) {
    ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
    
    for(i in 1:nrow(simil)){
      ad[which(simil[,i]>=z),i]<-1
      ad[which(simil[,i]<z),i]<-0
    }
    
    diag(ad)<-0
    gr=graph.adjacency(ad,mode="undirected",diag=FALSE)
    
    pvalue <- fit_power_law(degree(gr)) 
    
    if(pvalue$KS.p <= fit ){
      fit <- pvalue$KS.p
      value <- z                     
    }
  }
  
  Ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
  
  for(i in 1:nrow(simil)){
    Ad[which(simil[,i]>=value),i]<-1
    Ad[which(simil[,i]<value),i]<-0
  }
  
  colnames(Ad)<-rownames(Ad)<-rownames(simil)
  diag(Ad)<-0
  Gr=graph.adjacency(Ad,mode="undirected",add.colnames=NULL,diag=FALSE)
  print(paste("p-value =",fit,sep = " "))
  return(Gr) 
}

##################################################

Aarray <- getaffy(GSE = "GSETOTAL")
t <- c(rep(0,6),rep(1,12))
gene <- GeneSymbol("GPL570")
Adife <- difexprs(affy = Aarray,treatment = t,fdr = 0.45)
Anet <- CreateNet(difexp = Adife)
write.graph(Anet,file = "GSE29652NET.txt",format = "ncol")
