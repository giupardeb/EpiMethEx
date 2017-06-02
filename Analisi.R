start.time <- Sys.time()
library(data.table)
library(matrixStats)
#library(splitstackshape)
#library(doParallel)
#library(foreach)

#cl <- makeCluster(2)
dataIsLinearUser = T #da chiedere all'utente

setwd("/home/giuseppe/Scrivania/TesiRepoGithub/utilitiesTesi/")
source("functions.R")

dfPancan <-
  fread("dataset/pancan_normalized/genomicMatrix",
        header = T,
        sep = "\t")

#controllare cosa succede se viene eliminata la riga che contiene almeno un valore nella colonna
dfGpl <- na.omit(
  fread(
    "dataset/probe/GPL13534-11288.txt",
    sep = "\t",
    header = F ,
    skip = 37,
    na.strings = c("", "NA")
  )[-1, c(1, 22, 23, 24)]
)

dfMethylation <<-
  fread("dataset/Methylation450k/genomicMatrix", sep = "\t")
dfMethylation <- na.omit(dfMethylation)

# transpose all but the first column (name)
dfPancan2 <- as.data.frame(t(dfPancan[, -1]))
colnames(dfPancan2) <- dfPancan$sample

#rimuovo tutti quei geni che hanno valori uguali a zero
dfPancan2 <- dfPancan2[, which(!apply(dfPancan2 == 0, 2, all))]


indexTCGA <- data.frame(sapply(seq_along(dfPancan2), function(x) {
  row.names(dfPancan2[order(dfPancan2[x], decreasing = T), x, drop = F])
}))

colnames(indexTCGA) <- colnames(dfPancan2)


dfPancan2 <- apply(dfPancan2, 2, sort, decreasing = T)

quantili <-
  data.frame(colQuantiles(dfPancan2, probs = c(0.33, 0.66, 0.99)))
setnames(quantili, c("perc33", "perc66", "perc99"))
quantili <- as.data.frame(t(quantili))

divArray <- dim(dfPancan2[1, ])[2] / 3

meansUP <- apply(dfPancan2[1:158, ], 2, mean)
meansMID <- apply(dfPancan2[159:158 * 2, ], 2, mean)
meansDOWN <- apply(dfPancan2[317:473, ], 2, mean)
#inserisco alla fine della colonna di ogni gene la MEDIA DEL GRUPPO DOWN perché è sempre dispari
dfPancan2 <- rbind(dfPancan2, meansDOWN)
dfPancan2 <- rbind(dfPancan2, meansMID)
dfPancan2 <- rbind(dfPancan2, meansUP)
remove(meansUP, meansMID, meansDOWN)

dfPancan2 <- rbind(dfPancan2, quantili)
remove(quantili)

#Fold change delle varie combinazioni (up vs mid, up vs down, etc..)
ifelse(dataIsLinearUser, calcFCGeneLinear(), calcFCGeneLog())

dfTtest <- data.frame(
  matrix(NA, nrow = 6, ncol = dim(dfPancan2)[2]),
  row.names = c(
    "ttest_UPvsMID",
    "pvalue_UPvsMID",
    "ttest_UPvsDOWN",
    "pvalue_UPvsDOWN",
    "ttest_MIDvsDOWN",
    "pvalue_MIDvsDOWN"
  )
)
colnames(dfTtest) <- colnames(dfPancan2)
dimDFpancan <- dim(dfPancan)[2]

#registerDoParallel(cl)

for (k in 1:dimDFpancan) {
  #UPvsMID
  ttester(dfPancan2[1:158, k], dfPancan2[159:316, k], 1, 2)

  #UPvsDOWN
  ttester(dfPancan2[1:158, k], dfPancan2[317:474, k], 3, 4)

  #MIDvsDOWN
  ttester(dfPancan2[159:316, k], dfPancan2[317:474, k], 5, 6)
}

dfPancan2 <- rbind(dfPancan2, dfTtest)
remove(dfTtest)
#stopCluster(cl)

s <- strsplit(dfGpl$V22, split = ";")
GPL2 <-
  data.frame(V1 = rep(dfGpl$V1, sapply(s, length)), V22 = unlist(s))

s <- strsplit(dfGpl$V23, split = ";")
GPL3 <- data.frame(V2 = NA, V23 = unlist(s))

s <- strsplit(dfGpl$V24, split = ";")
GPL4 <- data.frame(V3 = NA, V24 = unlist(s))

righeCheTiServono1 <- cbind(GPL2, GPL3, GPL4)

righeCheTiServono1 <-
  Filter(function(x)
    ! (all(x == "")), righeCheTiServono1)

righeCheTiServono1 <-
  righeCheTiServono1[order(as.character(righeCheTiServono1$V22)),]

maxOccurence <-
  max(as.data.frame(table(unlist(
    righeCheTiServono1$V22
  )))$Freq)



##APPLY+FOREACH
DFfinale <-
  data.frame(matrix(NA, nrow = maxOccurence * 6, ncol = ncol(dfPancan2)))
DFtmp <- data.frame(matrix(NA, nrow = dim(indexTCGA)[1], ncol = 1))

colnames(DFfinale) <- colnames(indexTCGA)

z <- 1 #serve a continuare se ? lo stesso gene
genePrevious <- "" #evito la sovrapposizione sullo stesso gene

dimIndexTC <- dim(indexTCGA)[1]
colNamesTC <- colnames(indexTCGA)
colNamesDFM <- colnames(dfMethylation)


#registerDoParallel(cl)

tryCatch({
  
  apply(righeCheTiServono1, 1, function(i) {
    
    DFtmp <<- data.frame(matrix(NA, nrow = dimIndexTC, ncol = 1))
    
    ifelse(i["V22"] != genePrevious, z <<- 1, '')
    
    indexGeneColonnaDF2 <<- which(i["V22"] == colNamesTC)
    
    cgRigaDf3 <<-
      which(i["V1"] == as.character(dfMethylation$sample))
    
    #Abbiamo trovato che esistono CG che si riferiscono a geni che non esistono in dfpancan
    if (length(cgRigaDf3) != 0 && length(indexGeneColonnaDF2)!= 0) {
      
      for (j in 1:dimIndexTC) {
        
        TCindexDf3 <<-
          which(as.character(indexTCGA[j, indexGeneColonnaDF2]) == colNamesDFM)
        
        DFtmp[j, 1] <<-
          dfMethylation[cgRigaDf3, TCindexDf3, with = F][[1]]
        
      }
      
      meansDOWN <<- mean(DFtmp[317:473, ])
      DFtmp <<- rbind(DFtmp, meansDOWN)
      
      A <<- t.test(DFtmp[1:158, 1], DFtmp[159:316, 1],
                   var.equal = F)[c('statistic', 'p.value')]
      
      DFfinale[1:2, indexGeneColonnaDF2] <<-
        c(A$statistic, A$p.value)
      
      A <<- t.test(DFtmp[1:158, 1], DFtmp[317:474, 1],
                   var.equal = F)[c('statistic', 'p.value')]
      
      DFfinale[(z + 2):(z + 3), indexGeneColonnaDF2] <<-
        c(A$statistic, A$p.value)
      
      A <<- t.test(DFtmp[159:316, 1], DFtmp[317:474, 1],
                   var.equal = F)[c('statistic', 'p.value')]
      
      DFfinale[(z + 4):(z + 5), indexGeneColonnaDF2] <<-
        c(A$statistic, A$p.value)
      
    }
    z <<- z + 6
    
    genePrevious <<- i["V22"]
  })
  
},
error=function(cond) {
  message("Here's the original error message:")
  message(cond)
  message("\n questo è l'indice indexGeneColonnaDF2: ")
  message(indexGeneColonnaDF2)
  message("\n questo è l'indice cgRigaDf3: ")
  message(cgRigaDf3)
  message("\n questo è l'indice TCindexDf3: ")
  message(TCindexDf3)

}
)


remove(k, z, dfGpl, GPL2, GPL3, GPL4, A, DFtmp)
save.image(file = "AnalisiMultigenica1giugno.Rdata")

end.time <- Sys.time()
time.taken <- end.time - start.time
#print(end.time - start.time)
write(time.taken, file = "tempi")