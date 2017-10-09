library(data.table)
library(matrixStats)
library(plyr)
library(dplyr)
library(xlsx)
library(miscTools)
library(RcmdrMisc)
library(diptest)#usato solo per test

dataIsLinearUser = T #da chiedere all'utente

source("Functions.R")

dfPancan <-
  fread("dataset/pancan_normalized/genomicMatrix",
        header = T,
        sep = "\t")[c(11090,8410,13632), ]

dfGpl <-
  fread(
    "dataset/probe/GPL13534-11288.txt",
    sep = "\t",
    header = F ,
    skip = 37,
    na.strings = c("", "NA")
  )[-1, c(1, 15, 16, 22, 23, 24, 25, 26)]


# create a new column `island` with the two columns collapsed together
dfGpl$island <-
  apply(dfGpl[, c('V26', 'V25')] , 1 , paste , collapse = "_")

# remove the unnecessary columns
dfGpl <- dfGpl[,-c('V26', 'V25')]

dfMethylation <<-
  fread("dataset/Methylation450k/genomicMatrix", sep = "\t")

# transpose all but the first column (name)
dfPancan2 <- as.data.frame(t(dfPancan[,-1]))
colnames(dfPancan2) <- dfPancan$sample

#rimuovo tutti quei geni che hanno valori uguali a zero
dfPancan2 <- dfPancan2[, which(!apply(dfPancan2 == 0, 2, all))]


indexTCGA <- data.frame(sapply(seq_along(dfPancan2), function(x) {
  row.names(dfPancan2[order(dfPancan2[x], decreasing = T), x, drop = F])
}))

colnames(indexTCGA) <- colnames(dfPancan2)
#indexTCGA<-read.table("indexTCGA.csv")
#file.remove("indexTCGA.csv")

dfPancan2 <- apply(dfPancan2, 2, sort, decreasing = T)

quantili <-
  data.frame(colQuantiles(dfPancan2, probs = c(0.33, 0.66, 0.99)))
setnames(quantili, c("perc33", "perc66", "perc99"))
quantili <- as.data.frame(t(quantili))

dfPancan2 <- as.data.frame(dfPancan2)
#ottengo l'indice della colonna che contiene tutti i valori diversi tra loro, 
#evito l'errore nel momento in cui utilizzo bin.var
index<-which.max(apply(dfPancan2, 2, function(x)length(unique(x))))[[1]]
dfPancan2$variable <-
  with(dfPancan2,
       bin.var(
         dfPancan2[, index],
         bins = 3,
         method = 'proportions',
         labels = c('D', 'M', 'UP')
       ))
meanUP <-
  apply(dfPancan2[which(dfPancan2$variable %in% "UP"), -ncol(dfPancan2)], 2, mean)
meanMID <-
  apply(dfPancan2[which(dfPancan2$variable %in% "M"), -ncol(dfPancan2)], 2, mean)
meanDOWN <-
  apply(dfPancan2[which(dfPancan2$variable %in% "D"), -ncol(dfPancan2)], 2, mean)

#inserisco alla fine della colonna di ogni gene la MEDIA DEL GRUPPO DOWN perché è sempre dispari
dfPancan2 <- rbind(dfPancan2, meanDOWN)
dfPancan2 <- rbind(dfPancan2, meanMID)
dfPancan2 <- rbind(dfPancan2, meanUP)
#remove(meanUP, meanMID, meanDOWN)
rm(meanUP, meanMID, meanDOWN)

dfPancan2 <- rbind.fill(dfPancan2, quantili)

#remove(quantili)
rm(quantili)
gc()

#Fold change delle varie combinazioni (up vs mid, up vs down, etc..)

dfPancan2 <- calcFC(dfPancan2, dataIsLinearUser)

dfTtest <- data.frame(matrix(NA,
                             nrow = 6,
                             ncol = dim(dfPancan2[, -ncol(dfPancan2)])[2]))
colnames(dfTtest) <- colnames(dfPancan2[,-ncol(dfPancan2)])
dimDFpancan <- dim(dfPancan2[,-ncol(dfPancan2)])[2]

for (k in 1:dimDFpancan) {
  #UPvsMID
  dfTtest[1:2, k] <-
    ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "M"), k])
  
  #UPvsDOWN
  dfTtest[3:4, k] <-
    ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
  
  #MIDvsDOWN
  dfTtest[5:6, k] <-
    ttester(dfPancan2[which(dfPancan2$variable %in% "M"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
  gc()
  
}

prova <- data.frame(matrix(NA,
                           nrow = 6,
                           ncol = dim(dfPancan2[, -ncol(dfPancan2)])[2]))

ncoll<-dim(dfPancan2[, -ncol(dfPancan2)])[2]
prova[1:ncoll]<-lapply(dfPancan2[1:ncoll],function(){
  ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "M"), k])
  ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
  ttester(dfPancan2[which(dfPancan2$variable %in% "M"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
  
})


dfPancan2 <- rbind.fill(dfPancan2, dfTtest)
remove(dfTtest)
#dfPancan2 <- read.csv("dfPancan2.csv",check.names = F)
row.names(dfPancan2)[476] <- "meanUP"
row.names(dfPancan2)[475] <- "meanMedium"
row.names(dfPancan2)[474] <- "meanDown"
row.names(dfPancan2)[477] <- "perc33"
row.names(dfPancan2)[478] <- "perc66"
row.names(dfPancan2)[479] <- "perc99"
row.names(dfPancan2)[480] <- "fc_UPvsMID"
row.names(dfPancan2)[481] <- "fc_UPvsDOWN"
row.names(dfPancan2)[482] <- "fc_MIDvsDOWN"
row.names(dfPancan2)[483] <- "ttest_UPvsMID"
row.names(dfPancan2)[484] <- "pvalue_UPvsMID"
row.names(dfPancan2)[485] <- "ttest_UPvsDOWN"
row.names(dfPancan2)[486] <- "pvalue_UPvsDOWN"
row.names(dfPancan2)[487] <- "ttest_MIDvsDOWN"
row.names(dfPancan2)[488] <- "pvalue_MIDvsDOWN"


s <- strsplit(dfGpl$V22, split = ";")
s1 <- strsplit(dfGpl$V23, split = ";")
s2 <- strsplit(dfGpl$V24, split = ";")

dfAnnotations <-
  data.frame(
    cg = rep(dfGpl$V1, sapply(s, length)),
    #V1
    V15 = rep(dfGpl$V15, sapply(s, length)),
    V16 = rep(dfGpl$V16, sapply(s, length)),
    gene = unlist(s),
    #V22
    V23 = unlist(s1),
    posizione = unlist(s2),
    #V24
    island = rep(dfGpl$island, sapply(s, length))
    
  )
#remove(s, s1, s2)
rm(s, s1, s2, dfGpl)
gc()
#rimuovo i geni e i cg che non sono all'interno di dfpancan2
dfAnnotations <-
  subset(dfAnnotations,
         as.character(dfAnnotations$gene) %in% names(dfPancan2))

dfMethylation <-
  subset(dfMethylation,
         as.character(dfMethylation$sample) %in% dfAnnotations$cg)

maxOccurence <-
  max(as.data.frame(table(unlist(dfAnnotations$gene)))$Freq)


#Riordino per ogni gene i cg relativi
dfAnnotations <- dfAnnotations %>% arrange(gene, V16)

#adesso creo un dataframe dei CG per ogni gene con l'ordine indicato da indexTCGA

#dfCGunique conterrà la coppia (CG-posizione) univoca xk ci sono cg ripetuti in cui cambia la colonna V23

dfCGunique <-
  dfAnnotations[!duplicated(dfAnnotations[, c(1, 6)]), ]

dfMethylation$sample <- as.factor(dfMethylation$sample)
dfMethylation <- as.data.frame(dfMethylation)
rownames(dfMethylation) <- dfMethylation$sample
setkey(as.data.table(dfMethylation), sample)
gc()
tmp <- tapply(lapply(1:nrow(dfCGunique),
                     function (i)
                       (dfMethylation[as.vector(dfCGunique[i, 1]), as.vector(indexTCGA[, as.vector(dfCGunique[i, "gene"])])])),
              factor(dfCGunique[, "gene"]), function (x) {
                unname(unlist(x))
              })
max.rows <- max(sapply(tmp, length))
DFCGorder <-
  do.call(cbind, lapply(tmp, function(x) {
    length(x) <- max.rows
    return (x)
  }))

### <<<<<<<<<<<<<<<FASE I:
### CREAZIONE, PER OGNI GENE, TABELLE CONTENENTI CG APPARTENENTI AI VARI GRUPPI CON CALCOLO FC PVALUE>>>>>>>

valExprGene <-
  dfPancan2[c(480:482, 484, 486, 488),-ncol(dfPancan2)]
posizioni <- as.vector(unique(dfCGunique$posizione))
geni <- colnames(DFCGorder)
isole <- as.vector(unique(dfCGunique$island))
isole <-
  Filter(function(x)
    ! any(grepl("NA_NA", x)), isole) #remove "NA_NA"
mFinaleCGglobali <- data.frame(matrix())
mFinaleCGposition <- data.frame(matrix())
mFinaleCGisland <- data.frame(matrix())
mFinaleCGunificati <- data.frame(matrix())

stratification <- as.data.frame(rep("UP", each = 158))
b <- as.data.frame(rep("Medium", each = 157))
c <- as.data.frame(rep("Down", each = 158))
names(stratification) <- "stratification"
names(b) <- "stratification"
names(c) <- "stratification"
stratification <- rbind(stratification, b)
stratification <- rbind(stratification, c)

nome_colonne_position <- ""
nome_colonne_island <- ""
nome_colonne_cgUniti <- ""
nome_colonne_cgGlobali <- ""

rm(b, c, dfPancan, dfMethylation, indexTCGA)
gc()
lengthGeni <- length(geni)
lengthPosizioni <- length(posizioni)
lengthIsole <- length(isole)

for (i in 1:lengthGeni) {
  m <-
    data.frame(matrix(
      DFCGorder[, i],
      nrow = 473,
      ncol = length(DFCGorder[, i]) / 473
    ))

  colnames(m) <-
    as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 1])
  
  
  keep.cols <- names(m) %in% NA
  m <- m [!keep.cols]
  colnames(m) <- paste(names(m), geni[i], sep = "_")
  
  ###GC globali
  m1 <- m
  #conservo gli indici delle colonne NA per poi reinserirli alla fine
  columnNA <- which(sapply(m1, function(x)
    all(is.na(x))))
  
  #rimuovo colonne NA
  m1 <- m1[, colSums(is.na(m1)) != nrow(m1)]
  
  if (length(m1) != 0) {
    m1 <- cbind(m1, stratification)
    m1 <- Analisi(m1)
    m1$stratification <- NULL
  }
  
  if (length(columnNA) != 0) {
    for (j in 1:length(columnNA)) {
      m1 <- insertCol(as.matrix(m1), columnNA[[j]], v = NA)
    }
  }
  
  colnames(m1) <-
    paste(as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 1]),
          as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 6]),
          geni[i],
          sep = "_")
 
  #prima di eliminare le righe dei valori calcolare la correlazione 
  #tra valore gene e singolo CG
  
  m1<- rbind(m1,cor(dfPancan2[c(1:473),geni[i]],m1[c(1:473),]))
  m1<- rbind(m1,cor(dfPancan2[c(1:473),geni[i]],m1[c(1:473),],method = "spearman"))
  #tale controllo viene fatto per evitare che una colonna di NA di m1 venga eliminata
  #del tutto. Così facendo non ho problemi con eventuali cbind.
  if (dim(m1)[1] > 473) {
    m1 <-
      as.data.frame(m1[-c(1:473),])
  } else{
    m1 <-
      as.data.frame(m1[-c(1:464),])
  }
  
  #inserisco i valori del gene di riferimento
  m1 <-
    data.frame(sapply(m1, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
  
  mFinaleCGglobali <- cbind(mFinaleCGglobali, m1)
  rm(m1)
  
  ###CG secondo la posizione nel gene (TS, Body ecc)
  mFinaleCGposition <-
    cbind(
      mFinaleCGposition,
      Analisi2(
        lengthPosizioni,
        i,
        posizioni,
        "posizione",
        nome_colonne_position,
        T
      )
    )
  
  ###CG secondo le isole CpG
  mFinaleCGisland <-
    cbind(mFinaleCGisland,
          Analisi2(lengthIsole, i, isole, "island", nome_colonne_island, F))
  
  ###CG unificati
  m4 <- m
  m4 <- as.data.frame(m4[, colSums(is.na(m4)) != nrow(m4)])
  if (dim(m4)[2] != 0) {
    if (dim(m4)[2] != 1) {
      m4 <- stack(m4)
      m4 <- as.matrix(m4[, -2])
    }
    num_row_m4 <- nrow(m4)
    m4 <- cbind(m4, stratification)
    colnames(m4) <- c("value", "stratification")
    m4 <- Analisi(m4)
    m4 <- as.data.frame(m4[, -2])
    
    nome_colonne_cgUniti <-
      paste(nome_colonne_cgUniti, paste("CG", geni[i], sep = "_"), sep = ",")
    #eliminare le righe dei valori dei cg
    m4 <- as.data.frame(m4[-c(1:num_row_m4), ])
    m4 <-
      data.frame(sapply(m4, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
    
    m4<- rbind(m4,cor(dfPancan2[c(1:473),geni[i]],m4[c(1:473),]))
    m4<- rbind(m4,cor(dfPancan2[c(1:473),geni[i]],m4[c(1:473),],method = "spearman"))
    mFinaleCGunificati <- cbind(mFinaleCGunificati, m4)
  }
  rm(m4)
  gc()
}



mFinaleCGglobali <-
  mFinaleCGglobali[,-1]

mFinaleCGglobali<-setRowNames(mFinaleCGglobali)

a <- strsplit(colnames(mFinaleCGglobali), split = "_")
a <- unlist(lapply(a, '[[', 1))
mFinaleCGglobali <-
  rbind(mFinaleCGglobali, as.character(dfCGunique[which(dfCGunique$cg %in% as.array(a)), 7]))
row.names(mFinaleCGglobali)[nrow(mFinaleCGglobali)] <- "island"

mFinaleCGglobali <- t(mFinaleCGglobali)

write.xlsx(mFinaleCGglobali, "CG_Globali.xlsx", sheetName = "Sheet1")

mFinaleCGposition <-
  mFinaleCGposition[,-1]

nome_colonne_position <- strsplit(nome_colonne_position, ",")[[1]]
colnames(mFinaleCGposition) <- nome_colonne_position[-1]

mFinaleCGposition<-setRowNames(mFinaleCGposition)

mFinaleCGposition <- t(mFinaleCGposition)

write.xlsx(mFinaleCGposition, "CG_poszione_gene.xlsx", sheetName = "Sheet1")


mFinaleCGunificati <-
  mFinaleCGunificati[,-1]

nome_colonne_cgUniti <- strsplit(nome_colonne_cgUniti, ",")[[1]]
colnames(mFinaleCGunificati) <- nome_colonne_cgUniti[-1]

mFinaleCGunificati<-setRowNames(mFinaleCGunificati)

mFinaleCGunificati <- t(mFinaleCGunificati)

write.xlsx(mFinaleCGunificati, "CG_unificati.xlsx", sheetName = "Sheet1")


mFinaleCGisland <-
  mFinaleCGisland[,-1]

nome_colonne_island <- strsplit(nome_colonne_island, ",")[[1]]
colnames(mFinaleCGisland) <- nome_colonne_island[-1]

mFinaleCGisland<-setRowNames(mFinaleCGisland)

mFinaleCGisland <- t(mFinaleCGisland)

write.xlsx(mFinaleCGisland, "CG_isole_gene.xlsx", sheetName = "Sheet1")


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<FASE II: BOXPLOT!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###
#ciclo in cui per ogni gene(colonna di DFCGorder) costruisco una matrice ad hoc per creare i boxplot


#CG per CpG_Island
# numrow = 473 * 20 #numeroCG * numerototalirighe
# 
# df <- matrix(NA, nrow = numrow, ncol = 3)
# colnames(df) <- c('cg', 'island', 'value')
# islands <-
#   as.character(dfAnnotations[which(dfAnnotations$gene %in% "MMP9"), 7])
# cg <- colnames(m)
# m1 <- subset(m, select = -c(stratification))
# m1 <- t(m1)
# 
# i <- 1
# for (j in 1:20) {
#   #identifica cg e islands
#   for (k in 1:473) {
#     df[i, 1] <- cg[j]
#     df[i, 2] <- islands[j]
#     df[i, 3] <- m1[j, k]
#     i <- i + 1
#   }
# }
# df <- as.data.frame(df)
# df$value <- as.numeric(as.character(df$value))
# df <- cbind(df, stratification = m$stratification)
# df$f1f2 <- interaction(df$cg, df$stratification)
# 
# write.csv(df, "prova.csv")
# 
# #cg per posizione
# 
# a <-
#   melt(setDT(as.data.frame(mFinaleCGglobali)),
#        id = 'id',
#        measure.vars = list(c(2, 4, 6), c(3, 5, 7)))[, variable := c('fc_UPvsMID', 'fc_UPvsDOWN', 'fc_MIDvsDOWN')[variable]][order(id)]
