library(affy)
library(siggenes)
library(GEOquery)
library(vsn)
library(igraph)
library(acde)

alzGSM <- read.table("Alzheimer_Chips.txt",stringsAsFactors = F)
prkGSM <- read.table("Parkinson_Chips.txt",stringsAsFactors = F)
GetInfo(prkGSM,"GPL570")
scmGSM <- read.table("MultipleSclerosis_Chips.txt",stringsAsFactors = F)
GetInfo(scmGSM,"GPL570")

############# FUNCION ##############

GeneSymbol <- function(GPL, d = "."){
  # Va al directorio donde esta el archivo GPL
  setwd(d)
  # Extrae la informacion del archivo .soft
  gpl <- getGEO(filename = paste0(GPL,".soft"))
  # Crea una tabla con todos los datos en el GPL
  sym <- Table(gpl)
  # Crea un data.frame con las sondas asociadas al gen al que corresponden
  ta <- data.frame(sym$ID, gsub(" /// ","-",sym$`Gene Symbol`), stringsAsFactors = F)
  names(ta) <- c("probe","gene")
  return(ta)
}

####################################

############# FUNCION ##############

GetInfo <- function(GSE,GPL,dir="."){
  
  setwd(dir)
  
  sapply(as.vector(t(GSE)), getGEOSuppFiles)
  
  files <- dir(".")[grep("^GSE[0-9]",dir("."),ignore.case = T)]
  
  for(j in files){
    untar(paste0(j,"/",j,"_RAW.tar"), exdir = paste0(j,"/"))
  }
  getGEOfile(GPL, destdir = ".")
}

####################################

############# FUNCION ##############

getaffy <- function(GSE){
  raw <- read.table(file = paste0(GSE,"/","filelist.txt"),sep = "\t",
                    header = T,comment.char = "", stringsAsFactors = F)
  GSMs <- raw$Name[grep(".CEL",raw$Name,ignore.case = T)]
  affy <- ReadAffy(filenames = as.character(GSMs), compress = T,
                   celfile.path = GSE)
  return(affy)
}

###################################

########## FUNCION ################

difexprs <- function(affy,treatment,fdr,NormalizeMethod,SummaryMethod,DifferentialMethod){
  
  medianProbe <- function(gene,array){
    marray <- as.data.frame(exprs(array))
    names(marray) <- gsub(".CEL.gz","",names(marray),ignore.case = T)
    uni <- data.frame(gene$gene,marray)
    wowithw <- uni[grep(paste0("^","$"),uni$gene.gene,ignore.case = T,invert = T),]
    
    g <- data.frame()
    
    for(i in unique(na.omit(wowithw$gene.gene))){
      
      a <- as.data.frame(t(sapply(wowithw[grep(paste0("^",i,"$"),wowithw$gene.gene),
                                          2:ncol(wowithw)],median)),row.names = i)
      g <- rbind.data.frame(g,a)
    }
    return(g)
  }
  
  ProbeFilter <- function(array,gpl){
    # Extrae los datos de la clase eset o Affybatch
    eset <- exprs(array)
    # Obtiene el promedio de los valores normalizados
    rows <- rowMeans(eset,na.rm = T)
    # Crea un data.frame con el nombre de la sonda, el gen al que corresponde y un valor 0
    # que sera reemplazado posteriormente.
    da <- data.frame(gene,stringsAsFactors = F)
    # cambio de nombres del data.frame
    names(da) <- c("a","b")
    # crea un data.frame con el nombre de la sonda asociado con su valor promedio.
    probemean <- data.frame(names(rows),rows, stringsAsFactors = F)
    # cambio de nombres del data.frame
    names(probemean) <- c("a","b")
    # Union de los dos data.frames anteriormente creados a partir de los nombres de las sondas
    merge <- merge.data.frame(probemean, da, by.x = "a", by.y = "a")
    # crea un data.frame organizado y listo para ser filtrado.
    dat <- data.frame(merge$b.y,merge$b.x, row.names = merge$a,stringsAsFactors = F)
    # Elimina los datos que no estan asociados a ningun gen
    total <- dat[-c(grep(paste0("^","$"),dat[,1])),]
    # Crea una lista de genes unicos
    onlygenes <- unique(total[,1])
    # Crea un data.frame vacio que sera llenado con los datos finales
    result <- data.frame()
    
    # Filtrado de genes.
    for(x in onlygenes){
      # Busca todos los datos asociados a un gen
      genes  <- total[grep(x,total[,1],fixed = T),]
      # Toma el valor maximo.
      max <- genes[grep(max(total[grep(paste0("^",x,"$"),total[,1]),2]),genes[,2]),]
      # Llenado del data.frame vacio
      if(length(result) == 0){
        result <- rbind(max)
      }else{
        result <- rbind(result,max)
      }
    }
    # convirtiendo datos normalizados de matrix a data.frame
    eset2 <- as.data.frame(eset)
    # Obteniendo los datos finales
    final <- eset2[row.names(result),]
    # Reemplazando los nombres de sondas por genes
    row.names(final) <- result$merge.b.y
    
    return(final)
  }
  
  if(NormalizeMethod == "vsn"){
    vsn <- expresso(affy,pmcorrect.method = "pmonly", bg.correct = F,
                           normalize.method = "vsn", summary.method = "avgdiff")
    print("summarizing")
    if(SummaryMethod == "Max"){
      eset <- ProbeFilter(vsn,gene) 
    }else if(SummaryMethod == "Median"){
      eset <- medianProbe(gene,vsn)
    }
  }else if(NormalizeMethod == "rma"){
    rma <- rma(affy) 
    print("summarizing")
    if(SummaryMethod == "Max"){
      eset <- ProbeFilter(rma,gene) 
    }else if(SummaryMethod == "Median"){
      eset <- medianProbe(gene,rma)
    }
  }
  print("Differential analysis")
  
  if(DifferentialMethod == "sam"){
    samr <- sam(data = eset,cl = treatment,B=100,rand=100)
    tab <- as.data.frame(samr@mat.fdr)
    tab <- tab[tab$FDR >= fdr,]
    if(nrow(tab) == 0){stop("No differentially expressed genes found")}
    value <- tab[nrow(tab),]
    plot(samr,value$Delta)
    sum <- summary(samr,value$Delta,entrez=F)
    dife <- sum@row.sig.genes
    genes <- eset[dife,]
    print(paste0("Achieved FDR: ",value$FDR))
  }else if(DifferentialMethod == "acde"){
    acde <- stp(eset,t,R = 100, PER = T,alpha = fdr)
    plot(acde)
    print(paste0("Achieved FDR: ",acde$astar))
    print(paste0("delta value: ",acde$tstar))
    list <- data.frame(acde$gNames, acde$dgenes)
    diff <- list[list$acde.dgenes != "no-diff.",]
    genes <- eset[diff$acde.gNames,]
  }
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
    
    diag(ady)<-0
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
  
  fit <-0.05
  
  for (z in as.vector(C)) {
    ad <- matrix(0,ncol = nrow(simil),nrow = nrow(simil))
    
    for(i in 1:nrow(simil)){
      ad[which(simil[,i]>=z),i]<-1
      ad[which(simil[,i]<z),i]<-0
    }
    
    diag(ad)<-0
    gr=graph.adjacency(ad,mode="undirected",diag=FALSE)
    
    pvalue <- fit_power_law(degree(gr))
    
    print(paste("PCC:",z))
    print(paste("p-value:",pvalue$KS.p))
    
    if(pvalue$KS.p >= fit ){
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
  
  print("")
  print("Final Network")
  print("")
  print(paste("p-value =",fit,sep = " "))
  print(paste("threshold =",value,sep = " "))
  
  return(Gr) 
}

##################################################

Aarray <- getaffy(GSE = "GSE14711")
t <- c(rep(2,3),rep(1,8))
gene <- GeneSymbol("GPL570")
Adife <- difexprs(affy = Aarray,treatment = t,fdr = 0.05,NormalizeMethod = "vsn",
                  SummaryMethod = "Median",DifferentialMethod = "acde")
write.table(row.names(Adife),file = "GSE14711GeneList_vma.txt",quote = F)
