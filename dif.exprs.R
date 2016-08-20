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
    
    # Differential analysis using am method
    
    samr <- sam(data = eset,cl = treatment,B=100,rand=100)
    
    # Obtain the fdr to different tresholds
    
    tab <- as.data.frame(samr@mat.fdr)
    
    # Filter the fdr values from the expected value
    
    tab <- tab[tab$FDR >= fdr,]
    
    # Menssage to empty result
    
    if(nrow(tab) == 0){stop("No differentially expressed genes found")}
    
    # Obtain the threshold value
    
    value <- tab[nrow(tab),]
    
    # Plotting the result of differential analisys
    
    plot(samr,value$Delta)
    
    # Summarizing the results of differential analisys
    
    sum <- summary(samr,value$Delta,entrez=F)
    
    # Obtain the names of genes differentially expressed
    
    dife <- sum@row.sig.genes
    
    # Filter the expression matrix with the genes diferentially expressed
    
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
