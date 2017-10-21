start.time <- Sys.time()
library(data.table)
library(matrixStats)
library(plyr)
library(dplyr)
library(xlsx)
library(miscTools)
library(RcmdrMisc)
library(diptest)
library(psych)
library(doParallel)
library(foreach)

dataIsLinearUser = F
no_cores <- detectCores()/2

source("Functions.R")

dfPancan <-
  fread("dataset/pancan_normalized/genomicMatrix",
        header = T,
        sep = "\t")[1001:2001, ]

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
dfGpl <- dfGpl[, -c('V26', 'V25')]

dfMethylation <<-
  fread("dataset/Methylation450k/genomicMatrix", sep = "\t")

# transpose all but the first column (name)
dfPancan2 <- as.data.frame(t(dfPancan[, -1]))
colnames(dfPancan2) <- dfPancan$sample

#remove all genes that has all values equal a zero
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

dfPancan2 <- as.data.frame(dfPancan2)

index <-
  which.max(apply(dfPancan2, 2, function(x)
    length(unique(x))))[[1]]
dfPancan2$variable <-
  with(dfPancan2,
       bin.var(
         dfPancan2[, index],
         bins = 3,
         method = 'proportions',
         labels = c('D', 'M', 'UP')
       ))
meanUP <-
  apply(dfPancan2[which(dfPancan2$variable %in% "UP"),-ncol(dfPancan2)], 2, mean)
meanMID <-
  apply(dfPancan2[which(dfPancan2$variable %in% "M"),-ncol(dfPancan2)], 2, mean)
meanDOWN <-
  apply(dfPancan2[which(dfPancan2$variable %in% "D"),-ncol(dfPancan2)], 2, mean)

dfPancan2 <- rbind(dfPancan2, meanDOWN, meanMID, meanUP)

rm(meanUP, meanMID, meanDOWN)

dfPancan2 <- rbind.fill(dfPancan2, quantili)

rm(quantili)
gc()

#Calculate Fold change of the combinations (up vs mid, up vs down, etc..)

dfPancan2 <- calcFC(dfPancan2, dataIsLinearUser)

# dfTtest <- data.frame(matrix(NA,
#                              nrow = 6,
#                              ncol = dim(dfPancan2[, -ncol(dfPancan2)])[2]))
# colnames(dfTtest) <- colnames(dfPancan2[,-ncol(dfPancan2)])
dimDFpancan <- dim(dfPancan2[, -ncol(dfPancan2)])[2]
#
# for (k in 1:dimDFpancan) {
#   #UPvsMID
#   dfTtest[1:2, k] <-
#     ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "M"), k])
#
#   #UPvsDOWN
#   dfTtest[3:4, k] <-
#     ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
#
#   #MIDvsDOWN
#   dfTtest[5:6, k] <-
#     ttester(dfPancan2[which(dfPancan2$variable %in% "M"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
#   gc()
#
# }

# Create cluster with desired number of cores
cl <- makeCluster(no_cores)

# Register cluster
registerDoParallel(cl)

resultUPvsMID <-
  foreach(k = 1:dimDFpancan, .combine = rbind) %dopar% {
    #UPvsMID
    ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "M"), k])
  }

resultUPvsDOWN <-
  foreach(k = 1:dimDFpancan, .combine = rbind) %dopar% {
    #UPvsDOWN
    ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
    
  }

resultMIDvsDOWN <-
  foreach(k = 1:dimDFpancan, .combine = rbind) %dopar% {
    #MIDvsDOWN
    ttester(dfPancan2[which(dfPancan2$variable %in% "M"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k])
  }
stopCluster(cl)

resultUPvsMID <- as.data.frame(t(resultUPvsMID))
resultUPvsDOWN <- as.data.frame(t(resultUPvsDOWN))
resultMIDvsDOWN <- as.data.frame(t(resultMIDvsDOWN))
colnames(resultUPvsMID) <- colnames(dfPancan2[, -ncol(dfPancan2)])
colnames(resultUPvsDOWN) <- colnames(dfPancan2[, -ncol(dfPancan2)])
colnames(resultMIDvsDOWN) <-colnames(dfPancan2[, -ncol(dfPancan2)])

dfPancan2 <-
  rbind.fill(dfPancan2, resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
remove(resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
gc()

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

rm(s, s1, s2, dfGpl)
gc()
#remove genes and cg aren't into dfPancan2
#rimuovo i geni e i cg che non sono all'interno di dfpancan2
dfAnnotations <-
  subset(dfAnnotations,
         as.character(dfAnnotations$gene) %in% names(dfPancan2))

dfMethylation <-
  subset(dfMethylation,
         as.character(dfMethylation$sample) %in% dfAnnotations$cg)

maxOccurence <-
  max(as.data.frame(table(unlist(dfAnnotations$gene)))$Freq)

#reorder genes related CG
#Riordino per ogni gene i cg relativi
dfAnnotations <- dfAnnotations %>% arrange(gene, V16)

#adesso creo un dataframe dei CG per ogni gene con l'ordine indicato da indexTCGA

#dfCGunique conterrà la coppia (CG-posizione) univoca xk ci sono cg ripetuti in cui cambia la colonna V23

dfCGunique <-
  dfAnnotations[!duplicated(dfAnnotations[, c(1, 6)]),]

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
  dfPancan2[c(480:482, 484, 486, 488), -ncol(dfPancan2)]
posizioni <- as.vector(unique(dfCGunique$posizione))
geni <- colnames(DFCGorder)
isole <- as.vector(unique(dfCGunique$island))
isole <-
  Filter(function(x)
    ! any(grepl("NA_NA", x)), isole) #remove "NA_NA"

stratification <- as.data.frame(rep("UP", each = 158))
b <- as.data.frame(rep("Medium", each = 157))
c <- as.data.frame(rep("Down", each = 158))
names(stratification) <- "stratification"
names(b) <- "stratification"
names(c) <- "stratification"
stratification <- rbind(stratification, b)
stratification <- rbind(stratification, c)

rm(b, c, dfPancan, dfMethylation, indexTCGA)
gc()
lengthGeni <- length(geni)
lengthPosizioni <- length(posizioni)
lengthIsole <- length(isole)


# Create cluster with desired number of cores
cl <- makeCluster(no_cores)

# Register cluster
registerDoParallel(cl)

clusterCall(cl, function()
  library(plyr))
clusterCall(cl, function()
  library(miscTools))
clusterCall(cl, function()
  library(psych))

###GC globali
mFinaleCGglobali <-
  foreach(i = 1:lengthGeni, .combine = cbind) %dopar% {
#for (i in 1:lengthGeni){
  flag<-F
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
    }else{
      flag<-T
    }
    
    if (length(columnNA) != 0) {
      for (j in 1:length(columnNA)) {
        m1 <- insertCol(as.matrix(m1), columnNA[[j]], v = NA)
      }
    }
    
    if(flag){
      df <- as.data.frame(matrix(NA,nrow = 9,ncol = dim(m1)[2]))
      m1<- rbind(m1,df)
    }
    
    #prima di eliminare le righe dei valori calcolare la correlazione
    #tra valore gene e singolo CG
    
    c <- as.data.frame(dfPancan2[c(1:473), geni[i]])
    mTmp <- as.data.frame(m1[c(1:473), ])
    a <- corr.test(c, mTmp, adjust = "none")
    m1 <-
      rbind(m1, as.numeric(a$r), as.numeric(a$p)) #add correlation and p-value
    
    #tale controllo viene fatto per evitare che una colonna di NA di m1 venga eliminata
    #del tutto. Così facendo non ho problemi con eventuali cbind.
    if (dim(m1)[1] > 473) {
      m1 <-
        as.data.frame(m1[-c(1:473), ])
    } else{
      m1 <-
        as.data.frame(m1[-c(1:464), ])
    }
    
    #inserisco i valori del gene di riferimento
    
    m1 <-
      data.frame(sapply(m1, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
    
    
    colnames(m1) <-
      paste(as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 1]),
            as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 6]),
            geni[i],
            sep = "_")
    
    #mFinaleCGglobali <- cbind(mFinaleCGglobali, m1)
    m1
    #rm(m1)
  }

# mFinaleCGglobali <-
#   mFinaleCGglobali[, -1]

mFinaleCGglobali <- setRowNames(mFinaleCGglobali)

a <- strsplit(colnames(mFinaleCGglobali), split = "_")
a <- unlist(lapply(a, '[[', 1))
mFinaleCGglobali <-
  rbind(mFinaleCGglobali, as.character(dfCGunique[which(dfCGunique$cg %in% as.array(a)), 7]))
row.names(mFinaleCGglobali)[nrow(mFinaleCGglobali)] <- "island"

mFinaleCGglobali <- t(mFinaleCGglobali)

write.xlsx(mFinaleCGglobali, "CG_Globali.xlsx", sheetName = "Sheet1")
remove(mFinaleCGglobali,a)
gc()
###CG unificati

mFinaleCGunificati <- foreach(i = 1:lengthGeni, .combine = cbind) %dopar% {

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
    
    
    m4 <- m
    m4 <- as.data.frame(m4[, colSums(is.na(m4)) != nrow(m4)])
    num_CG <- length(m4)
    
    if (dim(m4)[2] != 0) {
      if (dim(m4)[2] != 1) {
        m4 <- stack(m4)
        m4 <- as.matrix(m4[,-2])
      }
      num_row_m4 <- nrow(m4)
      m4 <- cbind(m4, stratification)
      colnames(m4) <- c("value", "stratification")
      m4 <- Analisi(m4)
      m4 <- as.data.frame(m4[,-2])
      
      # nome_colonne_cgUniti <-paste("CG", geni[i], sep = "_")
      #   #paste(nome_colonne_cgUniti, paste("CG", geni[i], sep = "_"), sep = ",")
      names(m4) <- "value"
      
      dfTmp <-
        as.data.frame(rep(dfPancan2[c(1:473), geni[i]], num_CG))
      m4Tmp <- as.data.frame(m4[c(1:num_row_m4), ])
      a <- corr.test(dfTmp, m4Tmp, adjust = "none")
      m4 <-
        rbind(m4, as.numeric(a$r), as.numeric(a$p)) #add correlation and p-value
      
      #eliminare le righe dei valori dei cg
      m4 <- as.data.frame(m4[-c(1:num_row_m4),])
      m4 <-
        data.frame(sapply(m4, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
      colnames(m4) <- paste("CG", geni[i], sep = "_")
      m4
      #mFinaleCGunificati <- cbind(mFinaleCGunificati, m4)
    }else{
       m4<-as.data.frame(matrix(NA,nrow = 17,ncol = 1))
       colnames(m4) <- paste("CG", geni[i], sep = "_")
       m4
     }
  }

# mFinaleCGunificati <-
#   mFinaleCGunificati[, -1]

# nome_colonne_cgUniti <- strsplit(nome_colonne_cgUniti, ",")[[1]]
# colnames(mFinaleCGunificati) <- nome_colonne_cgUniti[-1]

mFinaleCGunificati <- setRowNames(mFinaleCGunificati)

mFinaleCGunificati <- t(mFinaleCGunificati)

write.xlsx(mFinaleCGunificati, "CG_unificati.xlsx", sheetName = "Sheet1")

###CG secondo la posizione nel gene (TS, Body ecc)
mFinaleCGposition <-
  foreach(i = 1:lengthGeni, .combine = cbind) %dopar% {
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
    
    Analisi2(lengthPosizioni, i, posizioni, "posizione")
    
  }

# mFinaleCGposition <-
#   mFinaleCGposition[, -1]

# nome_colonne_position <- strsplit(nome_colonne_position, ",")[[1]]
# colnames(mFinaleCGposition) <- nome_colonne_position[-1]

mFinaleCGposition <- setRowNames(mFinaleCGposition)

mFinaleCGposition <- t(mFinaleCGposition)

write.xlsx(mFinaleCGposition, "CG_poszione_gene.xlsx", sheetName = "Sheet1")

###CG secondo le isole CpG
mFinaleCGisland <-
  foreach(i = 1:lengthGeni, .combine = cbind) %dopar% {
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
    
    Analisi2(lengthIsole, i, isole, "island")
    
  }

# mFinaleCGisland <-
#   mFinaleCGisland[, -1]

# nome_colonne_island <- strsplit(nome_colonne_island, ",")[[1]]
# colnames(mFinaleCGisland) <- nome_colonne_island[-1]

mFinaleCGisland <- setRowNames(mFinaleCGisland)

mFinaleCGisland <- t(mFinaleCGisland)

write.xlsx(mFinaleCGisland, "CG_isole_gene.xlsx", sheetName = "Sheet1")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken







# for (i in 1:lengthGeni) {
#   tryCatch({
#     m <-
#       data.frame(matrix(
#         DFCGorder[, i],
#         nrow = 473,
#         ncol = length(DFCGorder[, i]) / 473
#       ))
#
#     colnames(m) <-
#       as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 1])
#
#
#     keep.cols <- names(m) %in% NA
#     m <- m [!keep.cols]
#     colnames(m) <- paste(names(m), geni[i], sep = "_")
#
#     ###GC globali
#     m1 <- m
#     #conservo gli indici delle colonne NA per poi reinserirli alla fine
#     columnNA <- which(sapply(m1, function(x)
#       all(is.na(x))))
#
#     #rimuovo colonne NA
#     m1 <- m1[, colSums(is.na(m1)) != nrow(m1)]
#
#     if (length(m1) != 0) {
#       m1 <- cbind(m1, stratification)
#       m1 <- Analisi(m1)
#       m1$stratification <- NULL
#     }
#
#     if (length(columnNA) != 0) {
#       for (j in 1:length(columnNA)) {
#         m1 <- insertCol(as.matrix(m1), columnNA[[j]], v = NA)
#       }
#     }
#
#     colnames(m1) <-
#       paste(as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 1]),
#             as.character(dfCGunique[which(dfCGunique$gene %in% geni[i]), 6]),
#             geni[i],
#             sep = "_")
#
#     #prima di eliminare le righe dei valori calcolare la correlazione
#     #tra valore gene e singolo CG
#     c <- as.data.frame(dfPancan2[c(1:473), geni[i]])
#     mTmp <- as.data.frame(m1[c(1:473), ])
#     a <- corr.test(c, mTmp, adjust = "none")
#     m1 <-
#       rbind(m1, as.numeric(a$r), as.numeric(a$p)) #add correlation and p-value
#
#     #tale controllo viene fatto per evitare che una colonna di NA di m1 venga eliminata
#     #del tutto. Così facendo non ho problemi con eventuali cbind.
#     if (dim(m1)[1] > 473) {
#       m1 <-
#         as.data.frame(m1[-c(1:473), ])
#     } else{
#       m1 <-
#         as.data.frame(m1[-c(1:464), ])
#     }
#
#     #inserisco i valori del gene di riferimento
#     m1 <-
#       data.frame(sapply(m1, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
#
#     mFinaleCGglobali <- cbind(mFinaleCGglobali, m1)
#     rm(m1)
#
#     ###CG secondo la posizione nel gene (TS, Body ecc)
#     mFinaleCGposition <-
#       cbind(
#         mFinaleCGposition,
#         Analisi2(
#           lengthPosizioni,
#           i,
#           posizioni,
#           "posizione",
#           nome_colonne_position,
#           T
#         )
#       )
#
#     ###CG secondo le isole CpG
#     mFinaleCGisland <-
#       cbind(
#         mFinaleCGisland,
#         Analisi2(lengthIsole, i, isole, "island", nome_colonne_island, F)
#       )
#
#     ###CG unificati
#     m4 <- m
#     m4 <- as.data.frame(m4[, colSums(is.na(m4)) != nrow(m4)])
#     num_CG <- length(m4)
#
#     if (dim(m4)[2] != 0) {
#       if (dim(m4)[2] != 1) {
#         m4 <- stack(m4)
#         m4 <- as.matrix(m4[,-2])
#       }
#       num_row_m4 <- nrow(m4)
#       m4 <- cbind(m4, stratification)
#       colnames(m4) <- c("value", "stratification")
#       m4 <- Analisi(m4)
#       m4 <- as.data.frame(m4[,-2])
#
#       nome_colonne_cgUniti <-
#         paste(nome_colonne_cgUniti, paste("CG", geni[i], sep = "_"), sep = ",")
#       names(m4) <- "value"
#
#       dfTmp <-
#         as.data.frame(rep(dfPancan2[c(1:473), geni[i]], num_CG))
#       m4Tmp <- as.data.frame(m4[c(1:num_row_m4), ])
#       a <- corr.test(dfTmp, m4Tmp, adjust = "none")
#       m4 <-
#         rbind(m4, as.numeric(a$r), as.numeric(a$p)) #add correlation and p-value
#
#       #eliminare le righe dei valori dei cg
#       m4 <- as.data.frame(m4[-c(1:num_row_m4),])
#       m4 <-
#         data.frame(sapply(m4, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
#
#
#       mFinaleCGunificati <- cbind(mFinaleCGunificati, m4)
#     }
#     rm(m4)
#     gc()
#
#   }, error = function(error_condition) {
#     print(error_condition)
#     print("valore di i: ")
#     print(i)
#   })
#
# }