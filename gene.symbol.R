# gene.symbol
# Juan David Henao Sanchez
# Bioinformatics and Systems Biology | Universidad Nacional de Colombia

gene.symbol <- function(GPL, d = "."){
  # Extrae la informacion del archivo .soft
  gpl <- getGEO(filename = paste0(GPL,".soft"))
  # Crea una tabla con todos los datos en el GPL
  sym <- Table(gpl)
  # Crea un data.frame con las sondas asociadas al gen al que corresponden
  ta <- data.frame(sym$ID, gsub(" /// ","-",sym$`Gene Symbol`), stringsAsFactors = F)
  names(ta) <- c("probe","gene")
  return(ta)
}