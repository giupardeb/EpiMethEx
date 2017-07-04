getValueDFPancan2 <- function(nameGene, indexRowPancan2) {
  #dobbiamo scrivere LA PRIMA PARTE DEL GENE IN JSON.
  
  rowname <- rownames(dfPancan2)
  
  outputJson[[numGene]]["name"] <- nameGene
  
  length <- length(indexRowPancan2)
  
  for (z in 1:length) {
    nameRow <- rowname[indexRowPancan2[z]]
    
    outputJson[[numGene]][nameRow] <-
      dfPancan2[indexRowPancan2[z], nameGene]
    
  }
  
  assign('outputJson', outputJson, envir = .GlobalEnv)
}

getValueDFFinale <- function(nameGene, indexRowDfFinale) {
  dim <- length(which(righeCheTiServono1$V22 == nameGene))
  cgPrevious <<- ""
  
  if (dim > 1)
    #possono esistere Geni che non hanno CG
  {
    for (j in 1:dim) {
      nameCG <-
        as.character(righeCheTiServono1[which(righeCheTiServono1$V22 == nameGene), ]$V1[j])
      
      if (nameCG != cgPrevious)
      {
        #PRENDO LA POSIZIONE E L'IDENTIFICATIVO DEL CG DEL GENE CONSIDERATO
        AAA <-
          righeCheTiServono1[which(righeCheTiServono1$V1 == nameCG &
                                     righeCheTiServono1$V22 == nameGene), 3:4]
        
        outputJson[[numGene]][["CG"]][[nameCG]] <-
          c(valori = DFfinale[indexRowDfFinale + 6 * (j - 1), nameGene],
            apply(AAA, 1, function(x)
              paste0(x[1], ',', x[2])))
        #identificativo = as.character(AAA$V23),
        #posizione = as.character(AAA$V24))
        
        cgPrevious <- nameCG
      }
      
    }
  }
  else {
    outputJson[[numGene]][["CG"]]["errore"] <-
      "Questo gene non ha i corrispettivi CG"
  }
  
  assign('outputJson', outputJson, envir = .GlobalEnv)
}