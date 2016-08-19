library(affy)
library(siggenes)
library(GEOquery)
library(vsn)
library(igraph)
library(acde)
library(minet)

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
    treatment[treatment == 0] <- 2
    acde <- stp(eset,treatment,R = 100, PER = T,alpha = fdr)
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

CreateNet <- function(difexp, method){
  
  if(method == "corelation"){
    simil <- abs(cor(t(difexp),use =  "pairwise.complete.obs"))
  }else if(method == "mutual information"){
    presimil <- build.mim(t(difexp), estimator = "mi.shrink", disc = "globalequalwidth")
    simil<-sqrt(1-exp(-2*presimil))
    simil[which(is.na(simil))]<-0
  }
  
  
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
    
    print(paste0(val,"%"))
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

Aarray <- getaffy(GSE = "PRKTOTAL")
t <- as.vector(t(read.table("PRKTOTAL/treatment.txt")))
gene <- GeneSymbol("GPL570")
Adife <- difexprs(affy = Aarray,treatment = t,fdr = 0.05,NormalizeMethod = "vsn",
                  SummaryMethod = "Median",DifferentialMethod = "sam")

write.table(Adife,"PRKMATRIX.txt",quote = F)


Adife <- read.table("GSE13732MATRIX.txt",header = T)

net1 <- CreateNet(difexp = Adife,method = "corelation")
net2 <- CreateNet(difexp = Adife,method = "mutual information")

write.graph(net1,"GSE13732corelation.net",format = "ncol")
write.graph(net2,"GSE13732mutual_information.net",format = "ncol")

###################################################
################ BOXPLOT ########################
##############################################

cof.var <- function(genes,study,type,treatment, complete = FALSE){
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
  
  Aarray <- getaffy(GSE = study)
  
  vsn <- expresso(Aarray,pmcorrect.method = "pmonly", bg.correct = F,
                  normalize.method = "vsn", summary.method = "avgdiff")
  
  data <- medianProbe(gene = genes, array = vsn)
  
  names(data) <- treatment
  
  if(complete == FALSE){
    tdata <- data[names(data) == type]
  }else{
    tdata = data
  }
  
  
  tdata$mean <- rowMeans(tdata, na.rm = T)
  CV <- function(x){sd(x,na.rm = T)/mean(x,na.rm = T)}
  tdata$cv <- apply(tdata[,1:(ncol(tdata)-1)],1,CV)
  
  return(tdata)
}

gpl <- GeneSymbol(GPL = "GPL570")

test <- cof.var(genes = gpl,study = "GSE16759", treatment = t, type = "0")

############################# ALZHEIMER #################################

t <- c(rep(1,4),rep(0,4))

GSE16759_T <- cof.var(genes = gpl,study = "GSE16759", treatment = t, type = "1")
GSE16759_C <- cof.var(genes = gpl,study = "GSE16759", treatment = t, type = "0")
GSE16759_X <- cof.var(genes = gpl,study = "GSE16759", treatment = t, complete = TRUE)

t1 <- c(rep(1,6),rep(0,3))

GSE18309_T <- cof.var(genes = gpl,study = "GSE18309", treatment = t1, type = "1")
GSE18309_C <- cof.var(genes = gpl,study = "GSE18309", treatment = t1, type = "0")
GSE18309_X <- cof.var(genes = gpl,study = "GSE18309", treatment = t1, complete = T)

t2 <- c(rep(0,8),rep(1,22))

GSE28146_T <- cof.var(genes = gpl,study = "GSE28146", treatment = t2, type = "1")
GSE28146_C <- cof.var(genes = gpl,study = "GSE28146", treatment = t2, type = "0")
GSE28146_X <- cof.var(genes = gpl,study = "GSE28146", treatment = t2, complete = T)

t3 <- c(0,1,0,1)

GSE28379_T <- cof.var(genes = gpl,study = "GSE28379", treatment = t3, type = "1")
GSE28379_C <- cof.var(genes = gpl,study = "GSE28379", treatment = t3, type = "0")
GSE28379_X <- cof.var(genes = gpl,study = "GSE28379", treatment = t3, complete = T)

t4 <- c(rep(0,6),rep(1,12))

GSE29652_T <- cof.var(genes = gpl,study = "GSE29652", treatment = t4, type = "1")
GSE29652_C <- cof.var(genes = gpl,study = "GSE29652", treatment = t4, type = "0")
GSE29652_X <- cof.var(genes = gpl,study = "GSE29652", treatment = t4, complete = T)

t5 <- rep(c(0,1),10)

GSE4757_T <- cof.var(genes = gpl,study = "GSE4757", treatment = t5, type = "1")
GSE4757_C <- cof.var(genes = gpl,study = "GSE4757", treatment = t5, type = "0")
GSE4757_X <- cof.var(genes = gpl,study = "GSE4757", treatment = t5, complete = T)

t6 <- c(rep(0,173),rep(1,80))

GSE48350_T <- cof.var(genes = gpl,study = "GSE48350", treatment = t6, type = "1")
GSE48350_C <- cof.var(genes = gpl,study = "GSE48350", treatment = t6, type = "0")
GSE48350_X <- cof.var(genes = gpl,study = "GSE48350", treatment = t6, complete = T)

t7 <- c(rep(0,74),rep(1,87))

GSE5281_T <- cof.var(genes = gpl,study = "GSE5281", treatment = t7, type = "1")
GSE5281_C <- cof.var(genes = gpl,study = "GSE5281", treatment = t7, type = "0")
GSE5281_X <- cof.var(genes = gpl,study = "GSE5281", treatment = t7, complete = T)

t8<- rep(c(0,1,1),2)

GSE6276_T <- cof.var(genes = gpl,study = "GSE6276", treatment = t8, type = "1")
GSE6276_C <- cof.var(genes = gpl,study = "GSE6276", treatment = t8, type = "0")
GSE6276_X <- cof.var(genes = gpl,study = "GSE6276", treatment = t8, complete = T)

t9 <- c(rep(1,4),rep(0,4))

GSE66333_T <- cof.var(genes = gpl,study = "GSE66333", treatment = t9, type = "1")
GSE66333_C <- cof.var(genes = gpl,study = "GSE66333", treatment = t9, type = "0")
GSE66333_X <- cof.var(genes = gpl,study = "GSE66333", treatment = t9, complete = T)

################################# Esclerosis ##########################################

w <- c(rep(1,39),rep(0,30),rep(1,34),rep(0,10))

GSE13732_T <- cof.var(genes = gpl,study = "GSE13732", treatment = w, type = "1")
GSE13732_C <- cof.var(genes = gpl,study = "GSE13732", treatment = w, type = "0")
GSE13732_X <- cof.var(genes = gpl,study = "GSE13732", treatment = w, complete = T)

w1 <- c(rep(0,15),rep(1,15))

GSE14386_T <- cof.var(genes = gpl,study = "GSE14386", treatment = w1, type = "1")
GSE14386_C <- cof.var(genes = gpl,study = "GSE14386", treatment = w1, type = "0")
GSE14386_X <- cof.var(genes = gpl,study = "GSE14386", treatment = w1, complete = T)

w2 <- c(0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,1)

GSE16461_T <- cof.var(genes = gpl,study = "GSE16461", treatment = w2, type = "1")
GSE16461_C <- cof.var(genes = gpl,study = "GSE16461", treatment = w2, type = "0")
GSE16461_X <- cof.var(genes = gpl,study = "GSE16461", treatment = w2, complete = T)

w3 <- c(rep(0,15),rep(1,14))

GSE21942_T <- cof.var(genes = gpl,study = "GSE21942", treatment = w3, type = "1")
GSE21942_C <- cof.var(genes = gpl,study = "GSE21942", treatment = w3, type = "0")
GSE21942_X <- cof.var(genes = gpl,study = "GSE21942", treatment = w3, complete = T)

w4 <- c(rep(0,4),rep(1,6))

GSE23205_T <- cof.var(genes = gpl,study = "GSE23205", treatment = w4, type = "1")
GSE23205_C <- cof.var(genes = gpl,study = "GSE23205", treatment = w4, type = "0")
GSE23205_X <- cof.var(genes = gpl,study = "GSE23205", treatment = w4, complete = T)

w5 <- c(rep(1,6),rep(0,4))

GSE26484_T <- cof.var(genes = gpl,study = "GSE26484", treatment = w5, type = "1")
GSE26484_C <- cof.var(genes = gpl,study = "GSE26484", treatment = w5, type = "0")
GSE26484_X <- cof.var(genes = gpl,study = "GSE26484", treatment = w5, complete = T)

w6 <- c(rep(0,8),rep(1,18))

GSE37750_T <- cof.var(genes = gpl,study = "GSE37750", treatment = w6, type = "1")
GSE37750_C <- cof.var(genes = gpl,study = "GSE37750", treatment = w6, type = "0")
GSE37750_X <- cof.var(genes = gpl,study = "GSE37750", treatment = w6, complete = T)

w7 <- c(0,0,rep(1,5))

GSE38010_T <- cof.var(genes = gpl,study = "GSE38010", treatment = w7, type = "1")
GSE38010_C <- cof.var(genes = gpl,study = "GSE38010", treatment = w7, type = "0")
GSE38010_X <- cof.var(genes = gpl,study = "GSE38010", treatment = w7, complete = T)

w8 <- c(rep(0,10),rep(1,10))

GSE43591_T <- cof.var(genes = gpl,study = "GSE43591", treatment = w8, type = "1")
GSE43591_C <- cof.var(genes = gpl,study = "GSE43591", treatment = w8, type = "0")
GSE43591_X <- cof.var(genes = gpl,study = "GSE43591", treatment = w8, complete = T)

w9 <- rep(c(1,0),8)

GSE52139_T <- cof.var(genes = gpl,study = "GSE52139", treatment = w9, type = "1")
GSE52139_C <- cof.var(genes = gpl,study = "GSE52139", treatment = w9, type = "0")
GSE52139_X <- cof.var(genes = gpl,study = "GSE52139", treatment = w9, complete = T)

w10 <- rep(c(0,1),6)

GSE53716_T <- cof.var(genes = gpl,study = "GSE53716", treatment = w10, type = "1")
GSE53716_C <- cof.var(genes = gpl,study = "GSE53716", treatment = w10, type = "0")
GSE53716_X <- cof.var(genes = gpl,study = "GSE53716", treatment = w10, complete = T)

w11 <- c(rep(0,7),rep(1,15))

GSE59085_T <- cof.var(genes = gpl,study = "GSE59085", treatment = w11, type = "1")
GSE59085_C <- cof.var(genes = gpl,study = "GSE59085", treatment = w11, type = "0")
GSE59085_X <- cof.var(genes = gpl,study = "GSE59085", treatment = w11, complete = T)
