library(data.table)
library(matrixStats)
library(plyr)
library(dplyr)
library(xlsx)
library(miscTools)
library(RcmdrMisc)

dataIsLinearUser = F #da chiedere all'utente

source("Functions.R")

dfPancan <-
  fread("dataset/pancan_normalized/genomicMatrix",
        header = T,
        sep = "\t")[c(8411,8412),]

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

dfPancan2$variable <-
  with(dfPancan2,
       bin.var(
         dfPancan2[c(1:nrow(dfPancan2)),1],
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

#inserisco alla fine della colonna di ogni gene la MEDIA DEL GRUPPO DOWN perché è sempre dispari
dfPancan2 <- rbind(dfPancan2, meanDOWN)
dfPancan2 <- rbind(dfPancan2, meanMID)
dfPancan2 <- rbind(dfPancan2, meanUP)
remove(meanUP, meanMID, meanDOWN)

dfPancan2 <- rbind.fill(dfPancan2, quantili)

remove(quantili)
#Fold change delle varie combinazioni (up vs mid, up vs down, etc..)

dfPancan2<-calcFC(dfPancan2,dataIsLinearUser)


dfTtest <- data.frame(matrix(NA,
                             nrow = 6,
                             ncol = dim(dfPancan2[,-ncol(dfPancan2)])[2]))

colnames(dfTtest) <- colnames(dfPancan2[, -ncol(dfPancan2)])
dimDFpancan <- dim(dfPancan2[, -ncol(dfPancan2)])[2]

for (k in 1:dimDFpancan) {
  #UPvsMID
  ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "M"), k], 1, 2)
  
  #UPvsDOWN
  ttester(dfPancan2[which(dfPancan2$variable %in% "UP"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k], 3, 4)
  
  #MIDvsDOWN
  ttester(dfPancan2[which(dfPancan2$variable %in% "M"), k], dfPancan2[which(dfPancan2$variable %in% "D"), k], 5, 6)
}

dfPancan2 <- rbind.fill(dfPancan2, dfTtest)
remove(dfTtest)
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

righeCheTiServono1 <-
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

#rimuovo i geni e i cg che non sono all'interno di dfpancan2
righeCheTiServono1 <-
  subset(righeCheTiServono1,
         as.character(righeCheTiServono1$gene) %in% names(dfPancan2))

dfMethylation <-
  subset(dfMethylation,
         as.character(dfMethylation$sample) %in% righeCheTiServono1$cg)

maxOccurence <-
  max(as.data.frame(table(unlist(
    righeCheTiServono1$gene
  )))$Freq)
remove(s, s1, s2)

#Riordino per ogni gene i cg relativi
righeCheTiServono1 <- righeCheTiServono1 %>% arrange(gene, V16)

#adesso creo un dataframe dei CG per ogni gene con l'ordine indicato da indexTCGA

#dfCGunique conterrà la coppia (CG-posizione) univoca xk ci sono cg ripetuti in cui cambia la colonna V23
dfCGunique <-
  righeCheTiServono1[!duplicated(righeCheTiServono1[, c(1, 6)]),]

dfMethylation$sample <- as.factor(dfMethylation$sample)
dfMethylation <- as.data.frame(dfMethylation)
rownames(dfMethylation) <- dfMethylation$sample
setkey(as.data.table(dfMethylation), sample)

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
mFinaleCGglobali <- data.frame(matrix(NA, ncol = 1))
mFinaleCGposition <- data.frame(matrix(NA, ncol = 1))
mFinaleCGisland <- data.frame(matrix(NA, ncol = 1))
mFinaleCGunificati <- data.frame(matrix(NA, ncol = 1))

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

for (i in 1:length(geni)) {
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
  m1 <- cbind(m1, stratification)
  m1 <- Analisi(m1)
  m1$stratification <- NULL
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
  m1 <-
    as.data.frame(m1[-c(1:473), ])
  
  m1 <-
    data.frame(sapply(m1, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
  
  mFinaleCGglobali <- cbind(mFinaleCGglobali, m1)
  
  
  ###CG secondo la posizione nel gene (TS, Body ecc)
  for (k in 1:length(posizioni)) {
    cg <-
      as.vector(dfCGunique[which(dfCGunique$gene %in% geni[i] &
                                   dfCGunique$posizione %in% posizioni[k]), 1])
    if (length(cg) != 0) {
      cg <- paste(cg, geni[i], sep = "_")
      m2 <- as.data.frame(m[, cg])
      
      #eliminare colonne tutte NA
      m2 <- as.data.frame(m2[, colSums(is.na(m2)) != nrow(m2)])
      if (dim(m2)[2] != 0) {
        if (dim(m2)[2] != 1) {
          m2 <- stack(m2)
          m2 <- as.matrix(m2[,-2])
        }
        num_row <- nrow(m2)
        m2 <- cbind(m2, stratification)
        colnames(m2) <- c("value", "stratification")
        
        #considerare se è necessario fare l'arrange..
        if (dim(m2)[2] != 2) {
          m2 <- arrange(m2, desc(stratification))
        }
        
        m2 <- Analisi(m2)
        m2 <- as.data.frame(m2[,-2])
        
        nome_colonne_position <-
          paste(nome_colonne_position,
                paste(posizioni[k], geni[i], sep = "_"),
                sep = ",")
        
        #eliminare le righe dei valori dei cg
        m2 <- as.data.frame(m2[-c(1:num_row),])
        m2 <-
          data.frame(sapply(m2, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
        mFinaleCGposition <- cbind(mFinaleCGposition, m2)
      }
    }
  }
  
  ###CG secondo le isole CpG
  
  for (j in 1:length(isole)) {
    cg <-
      as.vector(dfCGunique[which(dfCGunique$gene %in% geni[i] &
                                   dfCGunique$island %in% isole[j]), 1])
    if (length(cg) != 0) {
      cg <- paste(cg, geni[i], sep = "_")
      m3 <- as.data.frame(m[, cg])
      
      #eliminare colonne tutte NA
      m3 <- as.data.frame(m3[, colSums(is.na(m3)) != nrow(m3)])
      if (dim(m3)[2] != 0) {
        if (dim(m3)[2] != 1) {
          m3 <- stack(m3)
          m3 <- as.matrix(m3[,-2])
        }
        num_row <- nrow(m3)
        m3 <- cbind(m3, stratification)
        colnames(m3) <- c("value", "stratification")
        
        if (dim(m3)[2] != 2) {
          m3 <- arrange(m2, desc(stratification))
        }
        
        m3 <- Analisi(m3)
        m3 <- as.data.frame(m3[,-2])
        
        nome_colonne_island <-
          paste(nome_colonne_island,
                paste(isole[j], geni[i], sep = "_"),
                sep = ",")
        
        #eliminare le righe dei valori dei cg
        m3 <- as.data.frame(m3[-c(1:num_row),])
        m3 <-
          data.frame(sapply(m3, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
        
        mFinaleCGisland <- cbind(mFinaleCGisland, m3)
      }
    }
  }
  
  ###CG unificati
  m4 <- m
  m4 <- as.data.frame(m4[, colSums(is.na(m4)) != nrow(m4)])
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
    
    nome_colonne_cgUniti <-
      paste(nome_colonne_cgUniti, paste("CG", geni[i], sep = "_"), sep = ",")
    #eliminare le righe dei valori dei cg
    m4 <- as.data.frame(m4[-c(1:num_row_m4),])
    m4 <-
      data.frame(sapply(m4, c, unlist(valExprGene[, geni[i]])), row.names = NULL)
    mFinaleCGunificati <- cbind(mFinaleCGunificati, m4)
  }
  
}

mFinaleCGglobali <-
  mFinaleCGglobali[, -1]

row.names(mFinaleCGglobali)[1] <- "meanDown"
row.names(mFinaleCGglobali)[2] <- "meanMedium"
row.names(mFinaleCGglobali)[3] <- "meanUP"
row.names(mFinaleCGglobali)[4] <- "fc_UPvsMID"
row.names(mFinaleCGglobali)[5] <- "fc_UPvsDOWN"
row.names(mFinaleCGglobali)[6] <- "fc_MIDvsDOWN"
row.names(mFinaleCGglobali)[7] <- "pvalue_UPvsMID"
row.names(mFinaleCGglobali)[8] <- "pvalue_UPvsDOWN"
row.names(mFinaleCGglobali)[9] <- "pvalue_MIDvsDOWN"
row.names(mFinaleCGglobali)[10] <- "fc_UPvsMID(gene)"
row.names(mFinaleCGglobali)[11] <- "fc_UPvsDOWN(gene)"
row.names(mFinaleCGglobali)[12] <- "fc_MIDvsDOWN(gene)"
row.names(mFinaleCGglobali)[13] <- "pvalue_UPvsMID(gene)"
row.names(mFinaleCGglobali)[14] <- "pvalue_UPvsDOWN(gene)"
row.names(mFinaleCGglobali)[15] <- "pvalue_MIDvsDOWN(gene)"

a <- strsplit(colnames(mFinaleCGglobali), split = "_")
a <- unlist(lapply(a, '[[', 1))
mFinaleCGglobali <-
  rbind(mFinaleCGglobali, as.character(dfCGunique[which(dfCGunique$cg %in% as.array(a)), 7]))
row.names(mFinaleCGglobali)[16] <- "island"

mFinaleCGglobali <- t(mFinaleCGglobali)

write.xlsx(mFinaleCGglobali, "CG_Globali.xlsx", sheetName = "Sheet1")

mFinaleCGposition <-
  mFinaleCGposition[, -1]

nome_colonne_position <- strsplit(nome_colonne_position, ",")[[1]]
colnames(mFinaleCGposition) <- nome_colonne_position[-1]
row.names(mFinaleCGposition)[1] <- "meanDown"
row.names(mFinaleCGposition)[2] <- "meanMedium"
row.names(mFinaleCGposition)[3] <- "meanUP"
row.names(mFinaleCGposition)[4] <- "fc_UPvsMID"
row.names(mFinaleCGposition)[5] <- "fc_UPvsDOWN"
row.names(mFinaleCGposition)[6] <- "fc_MIDvsDOWN"
row.names(mFinaleCGposition)[7] <- "pvalue_UPvsMID"
row.names(mFinaleCGposition)[8] <- "pvalue_UPvsDOWN"
row.names(mFinaleCGposition)[9] <- "pvalue_MIDvsDOWN"
row.names(mFinaleCGposition)[10] <- "fc_UPvsMID(gene)"
row.names(mFinaleCGposition)[11] <- "fc_UPvsDOWN(gene)"
row.names(mFinaleCGposition)[12] <- "fc_MIDvsDOWN(gene)"
row.names(mFinaleCGposition)[13] <- "pvalue_UPvsMID(gene)"
row.names(mFinaleCGposition)[14] <- "pvalue_UPvsDOWN(gene)"
row.names(mFinaleCGposition)[15] <- "pvalue_MIDvsDOWN(gene)"
mFinaleCGposition <- t(mFinaleCGposition)

write.xlsx(mFinaleCGposition, "CG_poszione_gene.xlsx", sheetName = "Sheet1")


mFinaleCGunificati <-
  mFinaleCGunificati[, -1]

nome_colonne_cgUniti <- strsplit(nome_colonne_cgUniti, ",")[[1]]
colnames(mFinaleCGunificati) <- nome_colonne_cgUniti[-1]
row.names(mFinaleCGunificati)[1] <- "meanDown"
row.names(mFinaleCGunificati)[2] <- "meanMedium"
row.names(mFinaleCGunificati)[3] <- "meanUP"
row.names(mFinaleCGunificati)[4] <- "fc_UPvsMID"
row.names(mFinaleCGunificati)[5] <- "fc_UPvsDOWN"
row.names(mFinaleCGunificati)[6] <- "fc_MIDvsDOWN"
row.names(mFinaleCGunificati)[7] <- "pvalue_UPvsMID"
row.names(mFinaleCGunificati)[8] <- "pvalue_UPvsDOWN"
row.names(mFinaleCGunificati)[9] <- "pvalue_MIDvsDOWN"
row.names(mFinaleCGunificati)[10] <- "fc_UPvsMID(gene)"
row.names(mFinaleCGunificati)[11] <- "fc_UPvsDOWN(gene)"
row.names(mFinaleCGunificati)[12] <- "fc_MIDvsDOWN(gene)"
row.names(mFinaleCGunificati)[13] <- "pvalue_UPvsMID(gene)"
row.names(mFinaleCGunificati)[14] <- "pvalue_UPvsDOWN(gene)"
row.names(mFinaleCGunificati)[15] <- "pvalue_MIDvsDOWN(gene)"
mFinaleCGunificati <- t(mFinaleCGunificati)

write.xlsx(mFinaleCGunificati, "CG_unificati.xlsx", sheetName = "Sheet1")


mFinaleCGisland <-
  mFinaleCGisland[, -1]

nome_colonne_island <- strsplit(nome_colonne_island, ",")[[1]]
colnames(mFinaleCGisland) <- nome_colonne_island[-1]
row.names(mFinaleCGisland)[1] <- "meanDown"
row.names(mFinaleCGisland)[2] <- "meanMedium"
row.names(mFinaleCGisland)[3] <- "meanUP"
row.names(mFinaleCGisland)[4] <- "fc_UPvsMID"
row.names(mFinaleCGisland)[5] <- "fc_UPvsDOWN"
row.names(mFinaleCGisland)[6] <- "fc_MIDvsDOWN"
row.names(mFinaleCGisland)[7] <- "pvalue_UPvsMID"
row.names(mFinaleCGisland)[8] <- "pvalue_UPvsDOWN"
row.names(mFinaleCGisland)[9] <- "pvalue_MIDvsDOWN"
row.names(mFinaleCGisland)[10] <- "fc_UPvsMID(gene)"
row.names(mFinaleCGisland)[11] <- "fc_UPvsDOWN(gene)"
row.names(mFinaleCGisland)[12] <- "fc_MIDvsDOWN(gene)"
row.names(mFinaleCGisland)[13] <- "pvalue_UPvsMID(gene)"
row.names(mFinaleCGisland)[14] <- "pvalue_UPvsDOWN(gene)"
row.names(mFinaleCGisland)[15] <- "pvalue_MIDvsDOWN(gene)"
mFinaleCGisland <- t(mFinaleCGisland)

write.xlsx(mFinaleCGisland, "CG_isole_gene.xlsx", sheetName = "Sheet1")


### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<FASE II: BOXPLOT!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>###
#ciclo in cui per ogni gene(colonna di DFCGorder) costruisco una matrice ad hoc per creare i boxplot


#CG per CpG_Island
numrow = 473 * 20 #numeroCG * numerototalirighe

df <- matrix(NA, nrow = numrow, ncol = 3)
colnames(df) <- c('cg', 'island', 'value')
islands <-
  as.character(righeCheTiServono1[which(righeCheTiServono1$gene %in% "MMP9"), 7])
cg <- colnames(m)
m1 <- subset(m, select = -c(stratification))
m1 <- t(m1)

i <- 1
for (j in 1:20) {
  #identifica cg e islands
  for (k in 1:473) {
    df[i, 1] <- cg[j]
    df[i, 2] <- islands[j]
    df[i, 3] <- m1[j, k]
    i <- i + 1
  }
}
df <- as.data.frame(df)
df$value <- as.numeric(as.character(df$value))
df <- cbind(df, stratification = m$stratification)
df$f1f2 <- interaction(df$cg, df$stratification)

write.csv(df, "prova.csv")

#cg per posizione


#codice per saverio
Z <-
  righeCheTiServono1[!duplicated(righeCheTiServono1[, 1]), c(1, 4)]

Z <- Z[which(Z$gene %in% "TP53"), ]
AAA <- dfMethylation[which(dfMethylation$sample %in% Z$cg), ]
write.xlsx(AAA, "CG_TP53.xlsx", sheetName = "Sheet1")



a <-
  melt(setDT(as.data.frame(mFinaleCGglobali)),
       id = 'id',
       measure.vars = list(c(2, 4, 6), c(3, 5, 7)))[, variable := c('fc_UPvsMID', 'fc_UPvsDOWN', 'fc_MIDvsDOWN')[variable]][order(id)]
