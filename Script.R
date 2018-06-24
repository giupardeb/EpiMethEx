# [START] these istructions will download and will load a packages automatically. The goal is to load all packages, installed or not.
list.of.packages = c("data.table","matrixStats","plyr","dplyr","xlsx","miscTools","RcmdrMisc","diptest","psych","doParallel","foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)
# [END]

# [START] library for LIMMA
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
biocLite("gplots")
biocLite("pamr")
# [END]

# [START] third-partY library
source("Functions.R")
source("DownloadDataset.R")
# [END]

dataIsLinearUser = F
no_cores <- detectCores() / 2

#nameFolderDest <- "/home/giuseppe/dest"
#dir.create(file.path(nameFolderDest), showWarnings = FALSE)

nameFD <- "1-10" #dire al prof se possibile mettere questo come output
dir.create(file.path(nameFolderDest, nameFD))

# [START]: 
# params[1] = probe Dataset
# params[2] = Pancan Dataset
# params[3] = Methylation Dataset
# params[4] = Tumore Name
params <- datasetDownload("TCGA.SKCM.sampleMap/HumanMethylation450","TCGA.SKCM.sampleMap/HiSeqV2_PANCAN")
# [END]

#paths of the datasets
PathDatasetPancan <- file.path(getwd(),"datasets",params[4],params[2])
PathDatasetMethylation <- file.path(getwd(),"datasets",params[4],params[3])
PathDatasetProbe <- file.path(getwd(),"datasets",params[4],params[1])


dfPancan <- read.table(gzfile(PathDatasetPancan),header = T)[200:20, ]

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
dfPancan2 <- calcFC(dfPancan2, dataIsLinearUser)

dimDFpancan <- dim(dfPancan2[, -ncol(dfPancan2)])[2]

# Create cluster with desired number of cores
cl <- makeCluster(no_cores)

# Register cluster
registerDoParallel(cl)

#[START] calculate t-test for all genes
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
for(i in 475:489){
  row.names(dfPancan2)[i] <- rowNames[j]
  print(paste("valore i: ", i))
  print(paste("valore j: ", j))
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

#remove genes and cg aren't into dfPancan2 from dfAnnotations and dfMethylation
dfAnnotations <- subset(dfAnnotations, as.character(dfAnnotations$gene) %in% names(dfPancan2))

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

valExprGene <- dfPancan2[c(480:482, 484, 486, 488), -ncol(dfPancan2)]
positions <- as.vector(unique(dfCGunique$position))
genes <- colnames(DFCGorder)
islands <- as.vector(unique(dfCGunique$island))
islands <-
  Filter(function(x)
    ! any(grepl("NA_NA", x)), islands) #remove islands that have values "NA_NA"

# stratificationUp <- as.data.frame(rep("UP", each = 158))
# stratificationMedium <- as.data.frame(rep("Medium", each = 158))
# stratificationDown <- as.data.frame(rep("Down", each = 158))
# names(stratificationUp) <- "stratification"
# names(stratificationMedium) <- "stratification"
# names(stratificationDown) <- "stratification"
# stratification <- rbind(stratificationUp,stratificationMedium)
# stratification <- rbind(stratification, stratificationDown)
# remove(stratificationMedium, stratificationDown, dfPancan, dfMethylation, indexTCGA)
stratification <- as.data.frame(rep(c("UP","Medium","Down"), each = 158))
colnames(stratification) <- "stratification"


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

#[START] in questo ciclo for verrà creato un dataframe, in cui all'interno di esso ci saranno tutte le medie,le beta difference 
#        e i pvalue delle varie stratificazioni(UP,MED,DOWN), di tutti i CG in tutte le posizioni di tutti i geni.
mFinaleCGglobali <-
  foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
    
    flag <- F
    
    m <- data.frame(matrix(DFCGorder[, i],nrow = lengthDfpancan,ncol = length(DFCGorder[, i]) / lengthDfpancan ))
    colnames(m) <-as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
    keep.cols <- names(m) %in% NA
    m <- m [!keep.cols]
    colnames(m) <- paste(names(m), genes[i], sep = "_")
    
    m1 <- m
    columnNA <- which(sapply(m1, function(x) all(is.na(x))))
    
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
    
    c <- as.data.frame(dfPancan2[c(1:lengthDfpancan), genes[i]])
    mTmp <- as.data.frame(m1[c(1:lengthDfpancan), ])
    a <- corr.test(c, mTmp, adjust = "none")
    m1 <- rbind(m1, as.numeric(a$r), as.numeric(a$p))
    
    if (dim(m1)[1] > lengthDfpancan) {
      m1 <- as.data.frame(m1[-c(1:lengthDfpancan), ])
    } else{
      m1 <- as.data.frame(m1[-c(1:464), ])
    }
    
    m1 <- data.frame(sapply(m1, c, unlist(valExprGene[, genes[i]])), row.names = NULL)
    
    colnames(m1) <-
      paste(as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1]),
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 6]),
            genes[i],
            sep = "_")
    m1
  }
#[END]

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
    
    m <- data.frame(matrix( DFCGorder[, i], nrow = lengthDfpancan, ncol = length(DFCGorder[, i]) / lengthDfpancan ))
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
      
      dfTmp <- as.data.frame(rep(dfPancan2[c(1:lengthDfpancan), genes[i]], num_CG))
      m4Tmp <- as.data.frame(m4[c(1:num_row_m4), ])
      a <- corr.test(dfTmp, m4Tmp, adjust = "none")
      m4 <- rbind(m4, as.numeric(a$r), as.numeric(a$p)) #add correlation and p-value
      
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
    
    m <-data.frame(matrix(DFCGorder[, i],nrow = lengthDfpancan,ncol = length(DFCGorder[, i]) / lengthDfpancan ))
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
    
    m <- data.frame(matrix(DFCGorder[, i],nrow = lengthDfpancan, ncol = length(DFCGorder[, i]) / lengthDfpancan ))
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
