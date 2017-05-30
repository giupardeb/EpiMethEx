start.time <- Sys.time()

library(jsonlite)
FCuser <- 1.0
GeniSelezionati <- c()

source("functions.R")
#il which ritorna 1,2,3 ciò su dfpancan2 si traducono come: 1=480, 2=481, 3=483
DFfilterGene <- data.frame(dfPancan2[480:482,] >= FCuser)
DFfilterGene <-
  DFfilterGene[, which(!apply(DFfilterGene == F, 2, all))]

if (length(GeniSelezionati) != 0)
{
  #filtro il dataframe in base ai geni selezionati dall'utente
  DFfilterGene <-
    DFfilterGene[, GeniSelezionati[which(GeniSelezionati %in% names(DFfilterGene))]]
  
  GeniNonSelezionati <-
    GeniSelezionati[which(!GeniSelezionati %in% names(DFfilterGene))]
}


dimDF <- ncol(DFfilterGene)
outputJson <- list()

etichette <<-
  c(
    "ttest_UPvsMID",
    "pvalue_UPvsMID",
    "ttest_UPvsDOWN",
    "pvalue_UPvsDOWN",
    "ttest_MIDvsDOWN",
    "pvalue_MIDvsDOWN"
  )

for (i in 1:dimDF) {
  #aggiungo un gene nel JSON
  numGene <<- paste("gene", i, sep = "")
  outputJson[[numGene]][["CG"]] <- list()
  outputJson[[numGene]] <- list()
  
  namesGene <- colnames(DFfilterGene[, i, drop = F])
  
  if (DFfilterGene[1, i] == "TRUE" &&
      DFfilterGene[2, i] == "TRUE" &&
      DFfilterGene[3, i] == "TRUE") {
    #ELABORO TUTTO
    
    getValueDFPancan2(namesGene, c(477:488))
    getValueDFFinale(namesGene, c(1:6))
    
  }
  else if (DFfilterGene[1, i] == "TRUE" &&
           DFfilterGene[2, i] == "TRUE") {
    #ELABORO UPvsMID E UPvsDOWN
    
    getValueDFPancan2(namesGene, c(477:481, 483:486))
    getValueDFFinale(namesGene, c(1:4))
    
  }
  else if (DFfilterGene[1, i] == "TRUE" &&
           DFfilterGene[3, i] == "TRUE") {
    #ELABORO UPvsMID E MIDvsDOWN
    
    getValueDFPancan2(namesGene, c(477:480, 482:484, 487, 488))
    getValueDFFinale(namesGene, c(1:2, 5:6))
  }
  else if (DFfilterGene[2, i] == "TRUE" &&
           DFfilterGene[3, i] == "TRUE") {
    #ELABORO UPvsDOWN E MIDvsDOWN
    
    getValueDFPancan2(namesGene, c(477:479, 481, 482, 485:488))
    getValueDFFinale(namesGene, c(3:6))
  }
  else if (DFfilterGene[1, i] == "TRUE") {
    #ELABORO UPvsMID
    
    getValueDFPancan2(namesGene, c(477:480, 483, 484))
    getValueDFFinale(namesGene, c(1:2))
  }
  else if (DFfilterGene[2, i] == "TRUE") {
    #ELABORO UPvsDOWN
    
    getValueDFPancan2(namesGene, c(477:479, 481, 485, 486))
    getValueDFFinale(namesGene, c(3:4))
  }
  else if (DFfilterGene[3, i] == "TRUE") {
    #ELABORO MIDvsDOWN
    
    getValueDFPancan2(namesGene, c(477:479, 482, 487, 488))
    getValueDFFinale(namesGene, c(5:6))
  }
}

write(toJSON(outputJson, pretty = TRUE, auto_unbox = TRUE), file = "test2.json")
end.time <- Sys.time()
print(end.time - start.time)