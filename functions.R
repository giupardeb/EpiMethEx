calcFCLinear <- function(a,b) {
  
  #df <- data.frame()
  
  #meanUp <- dfPancan2[476,]
  #meanMid <- dfPancan2[475,]
  #meanDown <- dfPancan2[474,]
  
  fc <- setFC(a,b)
  #fc_UPvsDOWN =  setFC(meanUp,meanDown)
  #fc_MIDvsDOWN = setFC(meanMid,meanDown)
  
  # df <- rbind(df, fc_UPvsMID)
  # df <- rbind(df, fc_UPvsDOWN)
  # df <- rbind(df, fc_MIDvsDOWN)
  # rownames(df) <- c("fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN")
  # dfPancan2 <- rbind(dfPancan2, df)
  # remove(df)
  # assign('dfPancan2', dfPancan2, envir = .GlobalEnv)
  return(fc)
  
}

calcFCGeneLog <- function() {
  fc_UPvsMID = 2 ^ (abs(dfPancan2[475,] - dfPancan2[476,])) #meansMID / meansUP
  fc_UPvsDOWN = 2 ^ (abs(dfPancan2[474,] - dfPancan2[476,])) #meansDOWN / meansUP
  fc_MIDvsDOWN =  2 ^ (abs(dfPancan2[474,] - dfPancan2[475,])) #meansDOWN / meansMID
  
  df <- rbind(df, fc_UPvsMID)
  df <- rbind(df, fc_UPvsDOWN)
  df <- rbind(df, fc_MIDvsDOWN)
  rownames(df) <- c("fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN")
  dfPancan2 <- rbind(dfPancan2, df)
  remove(df)
  assign('dfPancan2', dfPancan2, envir = .GlobalEnv)
  
}

checkSign <- function(a,b) {
  return (sign(a)==sign(b))
}
setFC <- function(a,b){
  
  max <- checkMax(a,b)
  min <- checkMin(a,b)
    
  if(checkSign(a,b)){
    #segni concordi
    fc <- max / min
  }
  else {
    fc <- max - min
  }
  
  if(min(a,b) < max(a,b) )
    fc <- 1*fc
  else
    fc <- -1*fc
  
  return(fc)
}

checkMax <- function(a,b) {
  return(max(abs(a),abs(b)))
}
checkMin <- function(a,b) {
  return(min(abs(a),abs(b)))
}

ttester <- function(array1, array2, start, end) {
  if (sd(mapply('-', array1, array2, SIMPLIFY = T)) != 0) {
    A <- t.test(array1, array2,
                var.equal = F)[c('statistic', 'p.value')]
    dfTtest[start:end, k] <- c(A$statistic, A$p.value)
    
    assign('dfTtest', dfTtest, envir = .GlobalEnv)
  }
  
}


getValueDFPancan2 <- function(nameGene, indexRowPancan2) {
  #dobbiamo scrivere LA PRIMA PARTE DEL GENE IN JSON.
  
  rowname <- rownames(dfPancan2)
  
  outputJson[[numGene]]["name"] <- nameGene
  
  length <- length(indexRowPancan2)
  
  for (z in 1:length) {
    nameRow <- rowname[indexRowPancan2[z]]
    
    outputJson[[numGene]][nameRow] <-
      dfPancan2[indexRowPancan2[z], nameGene]
    
  }
  
  assign('outputJson', outputJson, envir = .GlobalEnv)
}

getValueDFFinale <- function(nameGene, indexRowDfFinale) {
  
  dim <- length(which(righeCheTiServono1$V22 == nameGene))
  cgPrevious <<- ""
  
  if (dim > 1)
    #possono esistere Geni che non hanno CG
  {
    for (j in 1:dim) {
      nameCG <-
        as.character(righeCheTiServono1[which(righeCheTiServono1$V22 == nameGene),]$V1[j])
      
      if (nameCG != cgPrevious)
      {
        #PRENDO LA POSIZIONE E L'IDENTIFICATIVO DEL CG DEL GENE CONSIDERATO
        AAA <-
          righeCheTiServono1[which(righeCheTiServono1$V1 == nameCG &
                                     righeCheTiServono1$V22 == nameGene), 3:4]
        
        outputJson[[numGene]][["CG"]][[nameCG]] <-
          c(
            valori = DFfinale[indexRowDfFinale + 6 * (j - 1), nameGene],
            apply(AAA, 1, function(x)paste0(x[1],',',x[2]))
            #identificativo = as.character(AAA$V23),
            #posizione = as.character(AAA$V24)
          )
        
        cgPrevious <- nameCG
      }
      
    }
  }
  else {
    outputJson[[numGene]][["CG"]]["errore"] <-
      "Questo gene non ha i corrispettivi CG"
  }
  
  assign('outputJson', outputJson, envir = .GlobalEnv)
}

ordinamento <- function(df,nameDF) {
  
  A <- as.data.frame(t(df))
  colnames(A) <- c("FC", "pValue")
  A <- cbind(A, V22 = rownames(A))
  
  if (nrow(A) > limit) {
    
    A <- A[order(A$FC, decreasing = T), ]
    A <- rbind(head(A, limit / 2), tail(A, limit / 2))
    
  }
  
  C <- merge(A, righeCheTiServono1, sort = F)
  
  colnames(C) <-
    c("Gene", "FC", "pValue", "CG", "Identificativo", "Posizione")
  
  #raggruppa FC GENE PVALUE, i cg sono all'interno di posizione
  Z<-C %>%
    group_by(Gene=C$Gene,FC=C$FC,pValue=C$pValue) %>%
    summarise(Posizione=paste(CG,Identificativo,Posizione,sep=";;",collapse=";;"))
  
  df <- Z[order(Z$FC, decreasing = T), ]
  
  assign(nameDF, df, envir = .GlobalEnv)
}