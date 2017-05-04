#https://rexplorations.wordpress.com/2015/08/11/confidence-intervals-and-t-tests-in-r/
#args <- commandArgs(TRUE) args dovr? contenere il path inviato da php

start.time <- Sys.time()
 library(RJSONIO)
 library(data.table)
 setwd("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/")
 #dfPancan <- read.delim("dataset/pancan_normalized/genomicMatrix", header = T, sep = "\t")
 #save.image("Work.Rdata")

load("Work.Rdata")

source("functions.R")

errori = 0
data = I
json1 <- list()
dataIsLinearUser = T #da chiedere all'utente
foldChangeUser = 0.05 #utilizzato per il foldchange, da chiedere all'utente


for(i in 1:dim(dfPancan)[1]){
#i=14516
  #controllo se la riga considerata è tutta nulla
  if(any(dfPancan[i,-1] == 0) == FALSE) {
    
    nameGene <- as.character(dfPancan[i,1])
    rowGene <- dfPancan[i,-1]
    rowGene <- sort(rowGene[1,],decreasing = T)
    divArray <- dim(rowGene[1,])[2]/3
    
    percentileGene <- quantile(rowGene[1,-1], c(0.33,0.66,0.99),na.rm = TRUE) #,na.rm = TRUE
    
    
    #questa media la faccio qui perché se l'array è pari devo aggiungere tale media all'arrayDown
    meanDOWN <- rowMeans(rowGene[1,(round(divArray)*2+1):(round(divArray)*3-1)])[[1]]
    
    if(divArray<round(divArray)){
      #approssimo per eccesso segue che aggiungo la media alla fine
      rowGene <- cbind(rowGene, mean_down=meanDOWN)
    }
    divArray <- round(divArray)
    
    #Faccio la media
    meanUP <- rowMeans(rowGene[1, 1:divArray])[[1]]
    meanMID <- rowMeans(rowGene[1,(divArray+1):(divArray*2)])[[1]]
    
    rowGene <- splitdf(rowGene,3)
    #Fold change delle varie combinazioni (up vs mid, up vs down, etc..)
    if(dataIsLinearUser){
      
      fc_UPvsMID = meanMID/meanUP
      fc_UPvsDOWN = meanDOWN/meanUP
      fc_MIDvsDOWN = meanDOWN/meanMID
      
    }
    else{
      
      fc_UPvsMID = 2^(abs(meanMID-meanUP))
      fc_UPvsDOWN = 2^(abs(meanDOWN-meanUP))
      fc_MIDvsDOWN =  2^(abs(meanDOWN-meanMID))
     
    }
    
    fc_Test = 0
    
    if(fc_UPvsMID > foldChangeUser){
      
      fc_Test = fc_Test + 1
      
    }
    
    if(fc_UPvsDOWN > foldChangeUser){
      
      fc_Test = fc_Test + 2
      
    }
    
    if(fc_MIDvsDOWN > foldChangeUser){
      
      fc_Test = fc_Test + 4
      
    }
    
    if(fc_Test == 0) {
      data <- list(
        list(
          errore = "Il gene non rispetta il foldchange desiderato"
        ))
      json1[[nameGene]] <- data
    }
    else{
      rowGeneUP <- as.data.frame(rowGene[1])
      rowGeneMID <- as.data.frame(rowGene[2])
      rowGeneDOWN <- as.data.frame(rowGene[3])
      json1[[nameGene]] <- metilazione()
    }
    
  }else {
    data <- list(
      list(
        errore = "Il gene ha tutti i valori settati a zero"
      ))
    json1[[nameGene]] <- data
  }
}

exportJson <- toJSON(json1)
write(exportJson, "result.json")
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)