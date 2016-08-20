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
    
    # Differential analysis using sam method
    
    samr <- sam(data = eset,cl = treatment,B=100,rand=100)
    
    # Obtain the fdr to different thresholds
    
    tab <- as.data.frame(samr@mat.fdr)
    
    # Filter the fdr values from the expected value
    
    tab <- tab[tab$FDR >= fdr,]
    
    # Message to empty result
    
    if(nrow(tab) == 0){stop("No differentially expressed genes found")}
    
    # Obtain the threshold value
    
    value <- tab[nrow(tab),]
    
    # Showing the result of differential analysis
    
    plot(samr,value$Delta)
    
    # Summarizing the results of differential analysis
    
    sum <- summary(samr,value$Delta,entrez=F)
    
    # Obtain the names of genes differentially expressed
    
    dife <- sum@row.sig.genes
    
    # Filter the expression matrix with the genes differentially expressed
    
    genes <- eset[dife,]
    
    # Showing the achieved fdr
    
    print(paste0("Achieved FDR: ",value$FDR))
    
  }else if(DifferentialMethod == "acde"){
    
    # Change the value of cases and controls
    
    treatment[treatment == 0] <- 2
    
    # Differential analysis using acde method
    
    acde <- stp(eset,treatment,R = 100, PER = T,alpha = fdr)
    
    # Showing the result of differential analysis
    
    plot(acde)
    
    # Showing the achieved fdr
    
    print(paste0("Achieved FDR: ",acde$astar))
    
    # Showing the threshold value
    
    print(paste0("delta value: ",acde$tstar))
    
    # Create a data.frame object with the result of differential analysis
    
    list <- data.frame(acde$gNames, acde$dgenes)
    
    # Obtain the name of genes differentially expressed
    
    diff <- list[list$acde.dgenes != "no-diff.",]
    
    # Filter the expression matrix with the genes differentially expressed
    
    genes <- eset[diff$acde.gNames,]
  }
 
   return(genes)
}
