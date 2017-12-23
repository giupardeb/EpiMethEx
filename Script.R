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
no_cores <- detectCores() / 2

nameFolderDest <- "/home/giuseppe/dest"
dir.create(file.path(nameFolderDest), showWarnings = FALSE)
nameFD <- "1-1000"
dir.create(file.path(nameFolderDest, nameFD))

PathDatasetPancan <- "dataset/pancan_normalized/genomicMatrix"
PathDatasetProbe <- "dataset/probe/GPL13534-11288.txt"
PathDatasetMethylation <- "dataset/Methylation450k/genomicMatrix"

source("Functions.R")

dfPancan <-fread(PathDatasetPancan,header = T,sep = "\t")[1:1000, ]

dfGpl <-fread(PathDatasetProbe,sep = "\t",header = F ,skip = 37,na.strings = c("", "NA"))[-1, c(1, 15, 16, 22, 23, 24, 25, 26)]

# create a new column `island` with the two columns collapsed together
dfGpl$island <-
  apply(dfGpl[, c('V26', 'V25')] , 1 , paste , collapse = "_")

# remove the unnecessary columns
dfGpl <- dfGpl[, -c('V26', 'V25')]

dfMethylation <-fread(PathDatasetMethylation, sep = "\t")

# transpose all but the first column (name)
dfPancan2 <- as.data.frame(t(dfPancan[, -1]))
colnames(dfPancan2) <- dfPancan$sample

#remove all genes that has all values equal a zero
dfPancan2 <- dfPancan2[, which(!apply(dfPancan2 == 0, 2, all))]

indexTCGA <- data.frame(sapply(seq_along(dfPancan2), function(x) {
  row.names(dfPancan2[order(dfPancan2[x], decreasing = T), x, drop = F])
}))

colnames(indexTCGA) <- colnames(dfPancan2)

dfPancan2 <- apply(dfPancan2, 2, sort, decreasing = T)

quantili <- data.frame(colQuantiles(dfPancan2, probs = c(0.33, 0.66, 0.99)))
setnames(quantili, c("perc33", "perc66", "perc99"))
quantili <- as.data.frame(t(quantili))

dfPancan2 <- as.data.frame(dfPancan2)

index <- which.max(apply(dfPancan2, 2, function(x) length(unique(x))))[[1]]

dfPancan2$variable <- with(dfPancan2,bin.var(dfPancan2[, index],bins = 3,method = 'proportions',labels = c('D', 'M', 'UP')))

meanUP <- apply(dfPancan2[which(dfPancan2$variable %in% "UP"),-ncol(dfPancan2)], 2, mean)
meanMID <- apply(dfPancan2[which(dfPancan2$variable %in% "M"),-ncol(dfPancan2)], 2, mean)
meanDOWN <- apply(dfPancan2[which(dfPancan2$variable %in% "D"),-ncol(dfPancan2)], 2, mean)

dfPancan2 <- rbind(dfPancan2, meanDOWN, meanMID, meanUP)

rm(meanUP, meanMID, meanDOWN)

dfPancan2 <- rbind.fill(dfPancan2, quantili)
rm(quantili)
gc()

#Calculate Fold change of the combinations (up vs mid, up vs down, etc..)
dfPancan2 <- calcFC(dfPancan2, dataIsLinearUser)

dimDFpancan <- dim(dfPancan2[, -ncol(dfPancan2)])[2]

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
colnames(resultMIDvsDOWN) <- colnames(dfPancan2[, -ncol(dfPancan2)])

dfPancan2 <- rbind.fill(dfPancan2, resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
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
dfAnnotations <- subset(dfAnnotations, as.character(dfAnnotations$gene) %in% names(dfPancan2))

dfMethylation <- subset(dfMethylation, as.character(dfMethylation$sample) %in% dfAnnotations$cg)

maxOccurence <- max(as.data.frame(table(unlist(dfAnnotations$gene)))$Freq)

#reorder genes related CG
dfAnnotations <- dfAnnotations %>% arrange(gene, V16)

#adesso creo un dataframe dei CG per ogni gene con l'ordine indicato da indexTCGA

#dfCGunique conterrà la coppia (CG-posizione) univoca xk ci sono cg ripetuti in cui cambia la colonna V23

dfCGunique <- dfAnnotations[!duplicated(dfAnnotations[, c(1, 6)]),]

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

valExprGene <- dfPancan2[c(480:482, 484, 486, 488), -ncol(dfPancan2)]
positions <- as.vector(unique(dfCGunique$posizione))
genes <- colnames(DFCGorder)
islands <- as.vector(unique(dfCGunique$island))
islands <-
  Filter(function(x)
    ! any(grepl("NA_NA", x)), islands) #remove islands that have values "NA_NA"

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
lengthGens <- length(genes)
lengthPositions <- length(positions)
lengthIslands <- length(islands)

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

###CG global
mFinaleCGglobali <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    flag <- F
    m <- data.frame(matrix(DFCGorder[, i],nrow = 473,ncol = length(DFCGorder[, i]) / 473 ))
    colnames(m) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    keep.cols <- names(m) %in% NA
    m <- m [!keep.cols]
    colnames(m) <- paste(names(m), genes[i], sep = "_")
    
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
    } else{
      flag <- T
    }
    
    if (length(columnNA) != 0) {
      for (j in 1:length(columnNA)) {
        m1 <- insertCol(as.matrix(m1), columnNA[[j]], v = NA)
      }
    }
    
    if (flag) {
      df <- as.data.frame(matrix(NA, nrow = 9, ncol = dim(m1)[2]))
      m1 <- rbind(m1, df)
    }
    
    c <- as.data.frame(dfPancan2[c(1:473), genes[i]])
    mTmp <- as.data.frame(m1[c(1:473), ])
    a <- corr.test(c, mTmp, adjust = "none")
    m1 <- rbind(m1, as.numeric(a$r), as.numeric(a$p))
    
    if (dim(m1)[1] > 473) {
      m1 <- as.data.frame(m1[-c(1:473), ])
    } else{
      m1 <- as.data.frame(m1[-c(1:464), ])
    }
    
    #inserisco i valori del gene di riferimento
    m1 <- data.frame(sapply(m1, c, unlist(valExprGene[, genes[i]])), row.names = NULL)
    
    colnames(m1) <-
      paste(as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1]),
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 6]),
            genes[i],
            sep = "_")
    m1
  }

mFinaleCGglobali <- setRowNames(mFinaleCGglobali)

a <- strsplit(colnames(mFinaleCGglobali), split = "_")
a <- unlist(lapply(a, '[[', 1))
mFinaleCGglobali <-rbind(mFinaleCGglobali, as.character(dfCGunique[which(dfCGunique$cg %in% as.array(a)), 7]))
row.names(mFinaleCGglobali)[nrow(mFinaleCGglobali)] <- "island"

mFinaleCGglobali <- t(mFinaleCGglobali)

write.xlsx(
  mFinaleCGglobali,
  paste(nameFolderDest, nameFD, "CG_data_individually.xlsx", sep = "/"),
  sheetName = "Sheet1"
)
remove(mFinaleCGglobali, a)
gc()

mFinaleCGunificati <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    m <- data.frame(matrix( DFCGorder[, i], nrow = 473, ncol = length(DFCGorder[, i]) / 473 ))
    colnames(m) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    keep.cols <- names(m) %in% NA
    m <- m [!keep.cols]
    colnames(m) <- paste(names(m), genes[i], sep = "_")
    
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
      
      names(m4) <- "value"
      
      dfTmp <- as.data.frame(rep(dfPancan2[c(1:473), genes[i]], num_CG))
      m4Tmp <- as.data.frame(m4[c(1:num_row_m4), ])
      a <- corr.test(dfTmp, m4Tmp, adjust = "none")
      m4 <- rbind(m4, as.numeric(a$r), as.numeric(a$p)) #add correlation and p-value
      
      #eliminare le righe dei valori dei cg
      m4 <- as.data.frame(m4[-c(1:num_row_m4),])
      m4 <- data.frame(sapply(m4, c, unlist(valExprGene[, genes[i]])), row.names = NULL)
      colnames(m4) <- paste("CG", genes[i], sep = "_")
      m4
      
    } else{
      m4 <- as.data.frame(matrix(NA, nrow = 17, ncol = 1))
      colnames(m4) <- paste("CG", genes[i], sep = "_")
      m4
    }
  }

mFinaleCGunificati <- setRowNames(mFinaleCGunificati)

mFinaleCGunificati <- t(mFinaleCGunificati)

write.xlsx(
  mFinaleCGunificati,
  paste(nameFolderDest, nameFD, "CG_of_a_genes.xlsx", sep = "/"),
  sheetName = "Sheet1"
)

###CG group by gene position
mFinaleCGposition <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    m <-data.frame(matrix(DFCGorder[, i],nrow = 473,ncol = length(DFCGorder[, i]) / 473 ))
    colnames(m) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    keep.cols <- names(m) %in% NA
    m <- m [!keep.cols]
    colnames(m) <- paste(names(m), genes[i], sep = "_")
    
    Analisi2(lengthPositions, i, positions, "posizione")
    
  }


mFinaleCGposition <- setRowNames(mFinaleCGposition)

mFinaleCGposition <- t(mFinaleCGposition)

write.xlsx(
  mFinaleCGposition,
  paste(nameFolderDest, nameFD, "CG_by_position.xlsx", sep = "/") ,
  sheetName = "Sheet1"
)

###CG group by island CpG
mFinaleCGisland <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    m <- data.frame(matrix(DFCGorder[, i],nrow = 473,ncol = length(DFCGorder[, i]) / 473 ))
    colnames(m) <- as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    keep.cols <- names(m) %in% NA
    m <- m [!keep.cols]
    colnames(m) <- paste(names(m), genes[i], sep = "_")
    
    Analisi2(lengthIslands, i, islands, "island")
    
  }

mFinaleCGisland <- setRowNames(mFinaleCGisland)

mFinaleCGisland <- t(mFinaleCGisland)

write.xlsx(
  mFinaleCGisland,
  paste(nameFolderDest, nameFD, "CG_by_islands.xlsx", sep = "/") ,
  sheetName = "Sheet1"
)
