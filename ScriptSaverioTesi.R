library(data.table)
library(matrixStats)
library(dplyr)
#cl <- makeCluster(2)
dataIsLinearUser = F #da chiedere all'utente

#setwd("/home/giuseppe/Scrivania/TesiRepoGithub/utilitiesTesi/")
source("functions.R")

dfPancan <-
  fread("dataset/pancan_normalized/genomicMatrix",
        header = T,
        sep = "\t")[c(10692, 8410, 17451),]

dfGpl <- na.omit(
  fread(
    "dataset/probe/GPL13534-11288.txt",
    sep = "\t",
    header = F ,
    skip = 37,
    na.strings = c("", "NA")
  )[-1, c(1, 15, 16, 22, 23, 24, 25, 26)]
)

# create a new column `x` with the three columns collapsed together
dfGpl$island <- apply( dfGpl[ , c('V26','V25') ] , 1 , paste , collapse = "_" )

# remove the unnecessary columns
dfGpl <- dfGpl[ , !( names( dfGpl ) %in% c('V25','V26') ) ]

dfMethylation <<-
  fread("dataset/Methylation450k/genomicMatrix", sep = "\t")
dfMethylation <- na.omit(dfMethylation)

# transpose all but the first column (name)
dfPancan2 <- as.data.frame(t(dfPancan[,-1]))
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

#divArray <- dim(dfPancan2[1, ])[2] / 3 #perché è qui? xD

meanUP <- apply(dfPancan2[1:158,], 2, mean)
meanMID <- apply(dfPancan2[159:316,], 2, mean)
meanDOWN <- apply(dfPancan2[317:473,], 2, mean)
#inserisco alla fine della colonna di ogni gene la MEDIA DEL GRUPPO DOWN perché è sempre dispari
dfPancan2 <- rbind(dfPancan2, meanDOWN)
dfPancan2 <- rbind(dfPancan2, meanMID)
dfPancan2 <- rbind(dfPancan2, meanUP)
remove(meanUP, meanMID, meanDOWN)

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
dimDFpancan <- dim(dfPancan2)[2]

for (k in 1:dimDFpancan) {
  #UPvsMID
  ttester(dfPancan2[1:158, k], dfPancan2[159:316, k], 1, 2)
  
  #UPvsDOWN
  ttester(dfPancan2[1:158, k], dfPancan2[317:474, k], 3, 4)
  
  #MIDvsDOWN
  ttester(dfPancan2[159:316, k], dfPancan2[317:474, k], 5, 6)
  k <- k + 1
}

dfPancan2 <- rbind(dfPancan2, dfTtest)
remove(dfTtest)

s <- strsplit(dfGpl$V22, split = ";")
GPL2 <-
  data.frame(
    V1 = rep(dfGpl$V1, sapply(s, length)),
    V22 = unlist(s),
    V15 = rep(dfGpl$V15, sapply(s, length)),
    V16 = rep(dfGpl$V16, sapply(s, length)),
    V25 = rep(dfGpl$V25, sapply(s, length))
  )

s <- strsplit(dfGpl$V23, split = ";")
GPL3 <- data.frame(V2 = NA, V23 = unlist(s))

s <- strsplit(dfGpl$V24, split = ";")
GPL4 <- data.frame(V3 = NA, V24 = unlist(s))

righeCheTiServono1 <- cbind(GPL2, GPL3, GPL4)

righeCheTiServono1 <-
  Filter(function(x)
    ! (all(x == "")), righeCheTiServono1)

righeCheTiServono1 <-
  righeCheTiServono1[order(as.character(righeCheTiServono1$V22)), ]

#rimuovo i geni e i cg che non sono all'interno di dfpancan2
righeCheTiServono1 <-
  subset(righeCheTiServono1,
         as.character(righeCheTiServono1[, 2]) %in% names(dfPancan2))

maxOccurence <-
  max(as.data.frame(table(unlist(
    righeCheTiServono1$V22
  )))$Freq)
remove(s)

#sarebbe bene chel'ordinamento dei CG per ogni gene avvenisse prima di fare l'analisi
righeCheTiServono1 <- righeCheTiServono1 %>% arrange(V16, V22)

##APPLY+FOREACH
DFfinale <<-
  data.frame(matrix(NA, nrow = maxOccurence * 6, ncol = ncol(dfPancan2)))
DFtmp <<- data.frame(matrix(NA, nrow = dim(indexTCGA)[1], ncol = 1))

colnames(DFfinale) <- colnames(indexTCGA)

z <<- 1 #serve a continuare se ? lo stesso gene
genePrevious <<- "" #evito la sovrapposizione sullo stesso gene

dimIndexTC <- dim(indexTCGA)[1]
colNamesTC <- colnames(indexTCGA)
colNamesDFM <- colnames(dfMethylation)

remove(dfGpl, GPL2, GPL3, GPL4, dfPancan)
#registerDoParallel(cl)

tryCatch({
  apply(righeCheTiServono1, 1, function(i) {
    DFtmp <<- matrix(NA, nrow = dimIndexTC, ncol = 1)
    
    ifelse(i["V22"] != genePrevious, z <<- 1, '')
    
    indexGeneColonna <<- which(i["V22"] == colNamesTC)
    
    cgRiga <<-
      which(i["V1"] == as.character(dfMethylation$sample))
    
    #esistono CG che si riferiscono a geni che non esistono in dfpancan
    if (length(cgRiga) != 0 && length(indexGeneColonna) != 0) {
      for (j in 1:dimIndexTC) {
        TCindex <<-
          which(as.character(indexTCGA[j, indexGeneColonna]) == colNamesDFM)
        
        DFtmp[j, 1] <<-
          dfMethylation[cgRiga, ..TCindex][[1]]
        
      }
      
      meanDOWN <<- mean(DFtmp[317:473, ])
      DFtmp <<- rbind(DFtmp, meanDOWN)
      
      #calcolo la media dei vari gruppi del CG
      meanUP <<- mean(DFtmp[1:158, ])
      meanMID <<- mean(DFtmp[159:316, ])
      meanDOWN <<- mean(DFtmp[317:474, ])
      
      if (dataIsLinearUser) {
        fcCG_UPvsMID <<- calcFCLinear(meanUP, meanMID)
        fcCG_UPvsDOWN <<- calcFCLinear(meanUP, meanDOWN)
        fcCG_MIDvsDOWN <<- calcFCLinear(meanMID, meanDOWN)
      } else{
        fcCG_UPvsMID <<- setFClog(meanUP, meanMID)
        fcCG_UPvsDOWN <<- setFClog(meanUP, meanDOWN)
        fcCG_MIDvsDOWN <<- setFClog(meanMID, meanDOWN)
      }
      
      
      A <<- t.test(DFtmp[1:158, 1], DFtmp[159:316, 1],
                   var.equal = F)['p.value']
      
      DFfinale[(z):(z + 1), indexGeneColonna] <<-
        c(fcCG_UPvsMID[[1]], A$p.value)
      
      A <<- t.test(DFtmp[1:158, 1], DFtmp[317:474, 1],
                   var.equal = F)['p.value']
      
      DFfinale[(z + 2):(z + 3), indexGeneColonna] <<-
        c(fcCG_UPvsDOWN[[1]], A$p.value)
      
      A <<- t.test(DFtmp[159:316, 1], DFtmp[317:474, 1],
                   var.equal = F)['p.value']
      
      DFfinale[(z + 4):(z + 5), indexGeneColonna] <<-
        c(fcCG_MIDvsDOWN[[1]], A$p.value)
      
    }
    z <<- z + 6
    
    genePrevious <<- i["V22"]
    
  })
  
},
error = function(cond) {
  message("Here's the original error message:")
  message(cond)
  message("\n questo è l'indice indexGeneColonna: ")
  message(indexGeneColonna)
  message("\n questo è l'indice cgRiga: ")
  message(cgRiga)
  message("\n questo è l'indice TCindex: ")
  message(TCindex)
  
})

#adesso creo un dataframe dei CG per ogni gene con l'ordine indicato da indexTCGA

#dfCGunique conterrà i CG univoci, in modo tale che DFCGorder contenga solo i CG univoci
dfCGunique <- righeCheTiServono1[!duplicated(righeCheTiServono1[, 1]), ]

dfMethylation$sample <- as.factor(dfMethylation$sample)
dfMethylation <- as.data.frame(dfMethylation)
rownames(dfMethylation) <- dfMethylation$sample
setkey(as.data.table(dfMethylation), sample)

tmp <- tapply(lapply(1:nrow(dfCGunique),
                     function (i)
                       (dfMethylation[as.vector(dfCGunique[i, 1]), as.vector(indexTCGA[, as.vector(dfCGunique[i, 2])])])),
              factor(dfCGunique[, 2]), function (x) {
                unname(unlist(x))
              })
max.rows <- max(sapply(tmp, length))
DFCGorder <-
  do.call(cbind, lapply(tmp, function(x) {
    length(x) <- max.rows
    return (x)
  }))

#ciclo in cui per ogni gene(colonna di DFCGorder) costruisco una matrice ad hoc per creare i boxplot

m <-
  data.frame(matrix(DFCGorder[, 1], nrow = 473, ncol = length(DFCGorder[, 1]) /
                      473))
colnames(m) <-
  as.character(dfCGunique[which(dfCGunique$V22 %in% "MMP9"), 1])

#cg con up mid down uniti
keep.cols <- names(m) %in% NA
m <- m [!keep.cols]
z <- head(rep(c("UP", "Medium", "Down"), each = 158), -1)
m$stratification <- z

#CG per posizione
numrow = 473 * 20 #numeroCG * numerototalirighe

df <- matrix(NA, nrow = numrow, ncol = 3)
colnames(df) <- c('cg', 'posizione', 'value')
posizioni <-
  as.character(righeCheTiServono1[which(righeCheTiServono1$V22 %in% "MMP9"), 7])
cg <- colnames(m)
m1 <- subset(m, select = -c(stratification))
m1 <- t(m1)

i <- 1
for (j in 1:20) {
  #identifica cg e posizioni
  for (k in 1:473) {
    df[i, 1] <- cg[j]
    df[i, 2] <- posizioni[j]
    df[i, 3] <- m1[j, k]
    i <- i + 1
  }
}
df <- as.data.frame(df)
df$value <- as.numeric(as.character(df$value))
df <- cbind(df, stratification = m$stratification)
df$f1f2 <- interaction(df$cg, df$stratification)

#boxplot CpG_Island_pos

s <-
  strsplit(as.character(righeCheTiServono1[which(righeCheTiServono1$V22 %in% "MMP9"), 5]), split = ":") #da inserire
namesColumns <-
  as.character(righeCheTiServono1[which(righeCheTiServono1$V22 %in% "MMP9"), 1])
GPL2 <-
  data.frame(V1 = rep(namesColumns, sapply(s, length)), V25 = unlist(s))#da inserire

nrows <- 473 * nrow(GPL2) #473=valori cg, 40 numero di righe di GPL2

df1 <- matrix(NA, nrow = nrows, ncol = 3)#da inserire
colnames(df1) <- c('cg', 'identificativo', 'value')#da inserire
posizioni <- as.character(GPL2$V25)#da inserire
cg <- namesColumns#da inserire
