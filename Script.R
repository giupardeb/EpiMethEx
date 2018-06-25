# [START] these istructions will download and will load a packages automatically. The goal is to load all packages, installed or not.
list.of.packages = c("data.table","matrixStats","plyr","dplyr","xlsx","miscTools","RcmdrMisc","diptest","psych","doParallel","foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)
# [END]

# # [START] library for LIMMA
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("limma")
# biocLite("gplots")
# biocLite("pamr")
# # [END]

# [START] third-partY library
source("Functions.R")
source("DownloadDataset.R")
# [END]

options(echo=TRUE) # To see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

#[START] input parameters
idMethylation <- args[1]
idPancan <- args[2]
minRangeGene <- args[3]
maxRangeGene <- args[4]
numCores <- as.numeric(args[5])
dataLinear <- as.logical(toupper(args[6]))
#[END]

#[START] static values, that identify the rows of DfPancan2
ROW_MEAN_DOWN <- 475
ROW_MEAN_MEDIUM <- 476
ROW_MEAN_UP <- 477
ROW_PERC_33 <- 478
ROW_PERC_66 <- 479
ROW_PERC_99 <- 480
ROW_FC_UPvsMID <- 481
ROW_FC_UPvsDOWN <- 482
ROW_FC_MIDvsDOWN <- 483
ROW_TTEST_UPvsMID <- 484
ROW_PVALUE_UPvsMID <- 485
ROW_TTEST_UPvsDOWN <- 486
ROW_PVALUE_UPvsDOWN <- 487
ROW_TTEST_MIDvsDOWN <- 488
ROW_PVALUE_MIDvsDOWN <- 489
#[END]

nameFolderDest <- file.path(getwd(),"outputs")
nameFD <- paste(minRangeGene,maxRangeGene,sep = "-")
dir.create(file.path(nameFolderDest, nameFD))

# [START]: 
# params[1] = probe Dataset
# params[2] = Pancan Dataset
# params[3] = Methylation Dataset
# params[4] = Tumore Name
params <- datasetDownload(idMethylation,idPancan)
# [END]

#paths of the datasets
PathDatasetPancan <- file.path(getwd(),"datasets",params[4],params[2])
PathDatasetMethylation <- file.path(getwd(),"datasets",params[4],params[3])
PathDatasetProbe <- file.path(getwd(),"datasets",params[4],params[1])

dfPancan <- read.table(gzfile(PathDatasetPancan),header = T)[as.numeric(minRangeGene):as.numeric(maxRangeGene), ]

#the first 37 rows of probset dataset, are comment. the following istruction read only columns of interest
dfGpl <-fread(PathDatasetProbe,sep = "\t",header = F ,skip = 37,na.strings = c("", "NA"))[-1, c(1, 15, 16, 22, 23, 24, 25, 26)]

# create a new column `island` with the two columns collapsed together
dfGpl$island <-
  apply(dfGpl[, c('V26', 'V25')] , 1 , paste , collapse = "_")

# remove the unnecessary columns
dfGpl <- dfGpl[, -c('V26', 'V25')]

dfMethylation <-read.table(gzfile(PathDatasetMethylation), sep = "\t", header = T)

# transpose all but the first column (name)
dfPancan2 <- as.data.frame(t(dfPancan[, -1]))
colnames(dfPancan2) <- dfPancan$sample

#remove all genes that has all values equal a zero
dfPancan2 <- dfPancan2[, which(!apply(dfPancan2 == 0, 2, all))]

#indexTCGA contains gene by gene the patient decreasing ordered
indexTCGA <- data.frame(sapply(seq_along(dfPancan2), function(x) {
  row.names(dfPancan2[order(dfPancan2[x], decreasing = T), x, drop = F])
}))

colnames(indexTCGA) <- colnames(dfPancan2)

dfPancan2 <- apply(dfPancan2, 2, sort, decreasing = T)

#calculate quantili of all genes into dfPancan2
quantili <- data.frame(colQuantiles(dfPancan2, probs = c(0.33, 0.66, 0.99)))
setnames(quantili, c("perc33", "perc66", "perc99"))
quantili <- as.data.frame(t(quantili))

dfPancan2 <- as.data.frame(dfPancan2)

# [START] these istruction divide into equal parts the dfPancan2
index <- which.max(apply(dfPancan2, 2, function(x) length(unique(x))))[[1]]
dfPancan2$variable <- with(dfPancan2,bin.var(dfPancan2[, index],bins = 3,method = 'proportions',labels = c('D', 'M', 'UP')))
# [END]
lengthDfpancan <- nrow(dfPancan2)
# [START] calculate the means of all values of genes that are UP, medium and down
meanUP <- apply(dfPancan2[which(dfPancan2$variable %in% "UP"),-ncol(dfPancan2)], 2, mean)
meanMID <- apply(dfPancan2[which(dfPancan2$variable %in% "M"),-ncol(dfPancan2)], 2, mean)
meanDOWN <- apply(dfPancan2[which(dfPancan2$variable %in% "D"),-ncol(dfPancan2)], 2, mean)
# [END]
dfPancan2 <- rbind(dfPancan2, meanDOWN, meanMID, meanUP)

remove(meanUP, meanMID, meanDOWN)

dfPancan2 <- rbind.fill(dfPancan2, quantili)
remove(quantili)
gc()

#Calculate Fold change of the combinations (up vs mid, up vs down, etc..)
dfPancan2 <- calcFC(dfPancan2, dataLinear)

dimDFpancan <- dim(dfPancan2[, -ncol(dfPancan2)])[2]

# Create cluster with desired number of cores
cl <- makeCluster(numCores)

# Register cluster
registerDoParallel(cl)

#[START] calculate t-test for all genes
resultUPvsMID <-
  foreach(k = 1:dimDFpancan, .combine = rbind) %dopar% {
    #UPvsMID
    calculateTtest(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "M"), k], T)
  }

resultUPvsDOWN <-
  foreach(k = 1:dimDFpancan, .combine = rbind) %dopar% {
    #UPvsDOWN
    calculateTtest(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k], T)
    
  }

resultMIDvsDOWN <-
  foreach(k = 1:dimDFpancan, .combine = rbind) %dopar% {
    #MIDvsDOWN
    calculateTtest(dfPancan2[which(dfPancan2$variable %in% "M"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k], T)
  }
stopCluster(cl)
#[END]

resultUPvsMID <- as.data.frame(t(resultUPvsMID))
resultUPvsDOWN <- as.data.frame(t(resultUPvsDOWN))
resultMIDvsDOWN <- as.data.frame(t(resultMIDvsDOWN))
colnames(resultUPvsMID) <- colnames(dfPancan2[, -ncol(dfPancan2)])
colnames(resultUPvsDOWN) <- colnames(dfPancan2[, -ncol(dfPancan2)])
colnames(resultMIDvsDOWN) <- colnames(dfPancan2[, -ncol(dfPancan2)])

dfPancan2 <- rbind.fill(dfPancan2, resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
remove(resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
gc()

rowNames <- c("meanDown","meanMedium","meanUP","perc33","perc66","perc99","fc_UPvsMID","fc_UPvsDOWN","fc_MIDvsDOWN",
              "ttest_UPvsMID","pvalue_UPvsMID","ttest_UPvsDOWN","pvalue_UPvsDOWN","ttest_MIDvsDOWN","pvalue_MIDvsDOWN")
j<-1
for(i in ROW_MEAN_DOWN:ROW_PVALUE_MIDvsDOWN){
  row.names(dfPancan2)[i] <- rowNames[j]
  j<-j+1
}

#[START] the dfAnnotations contains all cg with the corrispondent genes, islands, 
#        position of body and position into chromosoma
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
    position = unlist(s2),
    #V24
    island = rep(dfGpl$island, sapply(s, length))
  )
#[END]

remove(s, s1, s2, dfGpl)
gc()

#remove genes and cg aren't into dfPancan2 from dfAnnotations
dfAnnotations <- subset(dfAnnotations, as.character(dfAnnotations$gene) %in% names(dfPancan2))

#remove genes and cg aren't into dfPancan2 from dfMethylation
dfMethylation <- subset(dfMethylation, as.character(dfMethylation$sample) %in% dfAnnotations$cg)

maxOccurence <- max(as.data.frame(table(unlist(dfAnnotations$gene)))$Freq)

#reorder genes related CG
dfAnnotations <- dfAnnotations %>% arrange(gene, V16)
                         
#dfCGunique contain the unique couple (CG-position)
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

#DFCGorder is a dataframe that contains all CG values ordered by indexTCGA dataframe
DFCGorder <-
  do.call(cbind, lapply(tmp, function(x) {
    length(x) <- max.rows
    return (x)
  }))

valExprGene <- dfPancan2[c(ROW_FC_UPvsMID:ROW_FC_MIDvsDOWN, ROW_PVALUE_UPvsMID, ROW_PVALUE_UPvsDOWN, ROW_PVALUE_MIDvsDOWN), -ncol(dfPancan2)]
positions <- as.vector(unique(dfCGunique$position))
genes <- colnames(DFCGorder)
islands <- as.vector(unique(dfCGunique$island))
islands <-
  Filter(function(x)
    ! any(grepl("NA_NA", x)), islands) #remove islands that have values "NA_NA"

stratification <- as.data.frame(rep(c("UP","Medium","Down"), each = 158))
colnames(stratification) <- "stratification"

gc()
lengthGens <- length(genes)
lengthPositions <- length(positions)
lengthIslands <- length(islands)

# Create cluster with desired number of cores
cl <- makeCluster(numCores)

# Register cluster
registerDoParallel(cl)

clusterCall(cl, function()
  library(plyr))
clusterCall(cl, function()
  library(miscTools))
clusterCall(cl, function()
  library(psych))

#[START] in questo ciclo for verrà creato un dataframe, in cui all'interno di esso ci saranno tutte le medie,le beta difference 
#        e i pvalue delle varie stratificazioni(UP,MED,DOWN), di tutti i CG in tutte le posizioni di tutti i geni.
mFinaleCGglobali <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    flag <- F
    
    #[START] creo la matrice tempMatrix che conterrà i valori dei cg corrispondenti all'iesimo gene
    tempMatrix <- data.frame(matrix(DFCGorder[, i], nrow = lengthDfpancan, ncol = length(DFCGorder[, i]) / lengthDfpancan ))
    colnames(tempMatrix) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    #[END]
    
    #[START] elimino in tempMatrix le colonne che non sono associate a nessun cg dell'iesimo gene 
    keep.cols <- names(tempMatrix) %in% NA
    tempMatrix <- tempMatrix [!keep.cols]
    #[END]
    
    colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")
    
    #conservo su columnNA gli indici delle colonne che hanno tutti i valori NA
    columnNA <- which(sapply(tempMatrix, function(x) all(is.na(x))))
    
    tempMatrix <- tempMatrix[, colSums(is.na(tempMatrix)) != nrow(tempMatrix)]
    
    if (length(tempMatrix) != 0) {
      tempMatrix <- cbind(tempMatrix, stratification)
      tempMatrix <- Analysis(tempMatrix)
      tempMatrix$stratification <- NULL
    } else{
      flag <- T
    }
    
    #inserisco in tempMatrix le colonne che precedentemente ho eliminato nella stessa posizione
    if (length(columnNA) != 0) {
      for (j in 1:length(columnNA)) {
        tempMatrix <- insertCol(as.matrix(tempMatrix), columnNA[[j]], v = NA)
      }
    }
    
    if (flag) {
      df <- as.data.frame(matrix(NA, nrow = 9, ncol = dim(tempMatrix)[2]))
      tempMatrix <- rbind(tempMatrix, df)
    }
    
    #[START] calcolo la correlazione tra i dati di espressione dell'i-esimo gene con i dati di metilazione dei cg associati
    DataexpressionGeneTmp <- as.data.frame(dfPancan2[c(1:lengthDfpancan), genes[i]]) # prendo dati di espressione dei geni
    DataCG_Tmp <- as.data.frame(tempMatrix[c(1:lengthDfpancan), ]) #prendo i 474 valori di espressione dei cg
    resultCorrTest <- corr.test(DataexpressionGeneTmp, DataCG_Tmp , adjust = "none") #calcolo il test di correlazione tra i dati di espressione del gene e con i relativi CG
    tempMatrix <- rbind(tempMatrix, as.numeric(resultCorrTest$r), as.numeric(resultCorrTest$p))
    remove(DataexpressionGeneTmp,DataCG_Tmp,resultCorrTest)
    
    #[START] elimino da tempMatrix i valori di espressione dei cg, tenendo soltando le righe:
    #        "medianDown", "medianMedium", "medianUP", "bd_UPvsMID", "bd_UPvsDOWN", "bd_MIDvsDOWN", 
    #        "pvalue_UPvsMID","pvalue_UPvsDOWN", "pvalue_MIDvsDOWN", "pearson_correlation", "pvalue_pearson_correlation"
    
    if (dim(tempMatrix)[1] > lengthDfpancan) {
      tempMatrix <- as.data.frame(tempMatrix[-c(1:lengthDfpancan), ])
    } else{
      tempMatrix <- as.data.frame(tempMatrix[-c(1:464), ])
    }
    #[END]
    
    # insert into temp matrix "tempMatrix" the following values of the gene[i]:
    # "fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN", "pvalue_UPvsMID", "pvalue_UPvsDOWN", "pvalue_MIDvsDOWN"
    tempMatrix <- data.frame(sapply(tempMatrix, c, unlist(valExprGene[, genes[i]])), row.names = NULL)
    
    #[START] cambio il nome delle colonne associacndo il nome del cg con il corrispettivo posizione e gene.
    #        ex. cg27394127_TSS1500_BFAR
    colnames(tempMatrix) <-
      paste(as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1]),
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 6]),
            genes[i],
            sep = "_")
    #[END]
    
    #print the matrix1 into mFinaleCGglobali
    tempMatrix
  }
#[END]

mFinaleCGglobali <- setRowNames(mFinaleCGglobali)

#[START] estraggo i nomi dei cg e successivamente vado a recuperare tutte le isole associate a quei cg.
tmp <- strsplit(colnames(mFinaleCGglobali), split = "_")
tmp <- unlist(lapply(tmp, '[[', 1))

mFinaleCGglobali <-rbind(mFinaleCGglobali, as.character(dfCGunique[which(dfCGunique$cg %in% as.array(tmp)), 7]))
row.names(mFinaleCGglobali)[nrow(mFinaleCGglobali)] <- "island"
#[END]

mFinaleCGglobali <- t(mFinaleCGglobali)

write.xlsx(
  mFinaleCGglobali,
  file.path(nameFolderDest,nameFD,"CG_data_individually.xlsx"),
  #paste(nameFolderDest, nameFD, "CG_data_individually.xlsx", sep = "/"),
  sheetName = "Sheet1"
)
remove(mFinaleCGglobali, tmp)
gc()

mFinaleCGunificati <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    #[START] creo la matrice tempMatrix che conterrà i valori dei cg corrispondenti all'iesimo gene
    tempMatrix <- data.frame(matrix( DFCGorder[, i], nrow = lengthDfpancan, ncol = length(DFCGorder[, i]) / lengthDfpancan ))
    colnames(tempMatrix) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    #[END]
    
    #[START] elimino in tempMatrix le colonne che non sono associate a nessun cg dell'iesimo gene
    keep.cols <- names(tempMatrix) %in% NA
    tempMatrix <- tempMatrix [!keep.cols]
    #[END]
    
    colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")

    tempMatrix <- as.data.frame(tempMatrix[, colSums(is.na(tempMatrix)) != nrow(tempMatrix)])
    num_CG <- length(tempMatrix)
    
    if (dim(tempMatrix)[2] != 0) {
      if (dim(tempMatrix)[2] != 1) {
        tempMatrix <- stack(tempMatrix)
        tempMatrix <- as.matrix(tempMatrix[,-2])
      }
      num_row_m4 <- nrow(tempMatrix)
      tempMatrix <- cbind(tempMatrix, stratification)
      colnames(tempMatrix) <- c("value", "stratification")
      tempMatrix <- Analysis(tempMatrix)
      tempMatrix <- as.data.frame(tempMatrix[,-2])
      
      names(tempMatrix) <- "value"
      
      #[START] calcolo la correlazione tra i dati di espressione dell'i-esimo gene con i dati di metilazione dei cg associati
      dfTmp <- as.data.frame(rep(dfPancan2[c(1:lengthDfpancan), genes[i]], num_CG))
      m4Tmp <- as.data.frame(tempMatrix[c(1:num_row_m4), ])
      resultCorrTest <- corr.test(dfTmp, m4Tmp, adjust = "none")
      tempMatrix <- rbind(tempMatrix, as.numeric(resultCorrTest$r), as.numeric(resultCorrTest$p))
      #[END]
      remove(dfTmp,m4Tmp,resultCorrTest)
      
      tempMatrix <- as.data.frame(tempMatrix[-c(1:num_row_m4),])
      tempMatrix <- data.frame(sapply(tempMatrix, c, unlist(valExprGene[, genes[i]])), row.names = NULL)
      colnames(tempMatrix) <- paste("CG", genes[i], sep = "_")
      tempMatrix
      
    } else{
      tempMatrix <- as.data.frame(matrix(NA, nrow = 17, ncol = 1))
      colnames(tempMatrix) <- paste("CG", genes[i], sep = "_")
      tempMatrix
    }
  }

mFinaleCGunificati <- setRowNames(mFinaleCGunificati)

mFinaleCGunificati <- t(mFinaleCGunificati)

write.xlsx(
  mFinaleCGunificati,
  file.path(nameFolderDest,nameFD,"CG_of_a_genes.xlsx"),
  paste(nameFolderDest, nameFD, "CG_of_a_genes.xlsx", sep = "/"),
  sheetName = "Sheet1"
)

###CG group by gene position
mFinaleCGposition <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    #[START] creo la matrice tempMatrix che conterrà i valori dei cg corrispondenti all'iesimo gene
    tempMatrix <-data.frame(matrix(DFCGorder[, i],nrow = lengthDfpancan,ncol = length(DFCGorder[, i]) / lengthDfpancan ))
    colnames(tempMatrix) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    #[END]
    
    #[START] elimino in tempMatrix le colonne che non sono associate a nessun cg dell'iesimo gene
    keep.cols <- names(tempMatrix) %in% NA
    tempMatrix <- tempMatrix [!keep.cols]
    #[END]
    
    colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")
    
    AnalysisIslands_PositionsCG(lengthPositions, i, positions, "posizione")
    
  }

mFinaleCGposition <- setRowNames(mFinaleCGposition)

mFinaleCGposition <- t(mFinaleCGposition)

write.xlsx(
  mFinaleCGposition,
  file.path(nameFolderDest,nameFD,"CG_by_position.xlsx"),
  #paste(nameFolderDest, nameFD, "CG_by_position.xlsx", sep = "/") ,
  sheetName = "Sheet1"
)

###CG group by island CpG
mFinaleCGisland <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    #[START] creo la matrice tempMatrix che conterrà i valori dei cg corrispondenti all'iesimo gene
    tempMatrix <- data.frame(matrix(DFCGorder[, i],nrow = lengthDfpancan, ncol = length(DFCGorder[, i]) / lengthDfpancan ))
    colnames(tempMatrix) <- as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    #[END]
    
    #[START] elimino in tempMatrix le colonne che non sono associate a nessun cg dell'iesimo gene
    keep.cols <- names(tempMatrix) %in% NA
    tempMatrix <- tempMatrix [!keep.cols]
    #[END]
    
    colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")
    
    AnalysisIslands_PositionsCG(lengthIslands, i, islands, "island")
    
  }

mFinaleCGisland <- setRowNames(mFinaleCGisland)

mFinaleCGisland <- t(mFinaleCGisland)

write.xlsx(
  mFinaleCGisland,
  file.path(nameFolderDest,nameFD,"CG_by_islands.xlsx"),
  #paste(nameFolderDest, nameFD, "CG_by_islands.xlsx", sep = "/") ,
  sheetName = "Sheet1"
)
