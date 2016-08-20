# dif.exprs
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

dif.exprs <- function(affy,genes,treatment,fdr,NormalizeMethod,SummaryMethod,DifferentialMethod){
  
  if(NormalizeMethod == "vsn"){
    
    # Normalizing with vsn method
    
    vsn <- expresso(affy,pmcorrect.method = "pmonly", bg.correct = F,
                    normalize.method = "vsn", summary.method = "avgdiff")
    
    print("summarizing")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- .max.probe(vsn,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- .median.probe(gene,vsn)
    }
    
  }else if(NormalizeMethod == "rma"){
    
    # Normalizing using ram method
    
    rma <- rma(affy) 
    
    print("summarizing")
    
    if(SummaryMethod == "max"){
      
      # Summarizing using the high expression value
      
      eset <- ProbeFilter(rma,genes) 
      
    }else if(SummaryMethod == "median"){
      
      # Summarizing using the median expression value
      
      eset <- medianProbe(genes,rma)
    }
  }
  
  print("Differential analysis")
  
  if(DifferentialMethod == "sam"){
    
    # 
    
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
