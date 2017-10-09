calcBetaDifference <- function(matrix3){
  df <- data.frame()
  medianUP <- matrix3[nrow(matrix3), -ncol(matrix3)]
  medianMID <- matrix3[nrow(matrix3)-1, -ncol(matrix3)]
  medianDOWN <- matrix3[nrow(matrix3)-2, -ncol(matrix3)]
  
  bd_UPvsMID <- medianUP - medianMID
  bd_UPvsDOWN <- medianUP - medianDOWN
  bd_MIDvsDOWN <- medianMID - medianDOWN
  
  df <- rbind(df, bd_UPvsMID,bd_UPvsDOWN,bd_MIDvsDOWN)
  
  if(dim(matrix3)[2] == 2){
    colnames(df)<-"value"
  }else{
    colnames(df) <- colnames(matrix3[,-ncol(matrix3)])
  }
  
  matrix3 <- rbind.fill(matrix3,df)
  remove(df)
  return(matrix3)
}

calcFC <- function(matrix2,dataIsLinear) {
  df <- data.frame()
  
  meanUp <- matrix2[nrow(matrix2), -ncol(matrix2)]
  meanMid <- matrix2[nrow(matrix2)-1, -ncol(matrix2)]
  meanDown <- matrix2[nrow(matrix2)-2, -ncol(matrix2)]
  
  if(dataIsLinear){
    fc_UPvsMID <- setFC(meanUp, meanMid)
    fc_UPvsDOWN <-  setFC(meanUp, meanDown)
    fc_MIDvsDOWN <- setFC(meanMid, meanDown)
  }else{
    fc_UPvsMID <- setFClog(meanUp,meanMid) 
    fc_UPvsDOWN <- setFClog(meanUp,meanDown)
    fc_MIDvsDOWN <-  setFClog(meanMid,meanDown)
  }
  
  df <- rbind(df, fc_UPvsMID,fc_UPvsDOWN,fc_MIDvsDOWN)
  
  if(dim(matrix2)[2] == 2){
    colnames(df)<-"value"
  }else{
    colnames(df) <- colnames(matrix2[,-ncol(matrix2)])
  }
  
  matrix2 <- rbind.fill(matrix2,df)
  remove(df)
  return(matrix2)
}


checkSign <- function(a, b) {
  return (sign(a) == sign(b))
}

setFClog <- function(meanFirstGroup, meanSecondGroup){
  
  fold <- 2 ^ (abs(meanFirstGroup - meanSecondGroup))
  
  fc <- lapply(seq_along(fold), function(i) {
    if (meanSecondGroup[[i]] < meanFirstGroup[[i]]) {
      fold[[i]]
    }
    else{
      -fold[[i]]
    }
  })
  
  return(fc)
  
}
setFC <- function(meanFirstGroup, meanSecondGroup) {
  
  maxabs <- mapply(max, abs(meanFirstGroup), abs(meanSecondGroup))
  minabs <- mapply(min, abs(meanFirstGroup), abs(meanSecondGroup))
  max <- mapply(max, meanFirstGroup, meanSecondGroup)
  min <- mapply(min, meanFirstGroup, meanSecondGroup)
  Y <- checkSign(meanFirstGroup, meanSecondGroup)
  
  fc <- lapply(seq_along(Y), function(i) {
    if (Y[[i]] == T) {
      maxabs[[i]] / minabs[[i]]
    } else{
      max[[i]] - min[[i]]
    }
  })
  
  fold <- lapply(seq_along(fc), function(i) {
    if (meanSecondGroup[[i]] < meanFirstGroup[[i]]) {
      fc[[i]]
    }
    else{
      -fc[[i]]
    }
  })
  
  return(fold)
}

#Utilizzato nell'analisi dei Geni
ttester <- function(array1, array2) {
  
  if (sd(mapply('-', array1, array2, SIMPLIFY = T),na.rm = T) != 0) {
    
    A <- t.test(array1, array2,
                var.equal = F)[c('statistic', 'p.value')]
    return(c(A$statistic, A$p.value))
  }else{
    return(c(NA,NA))
  }
}

#Utilizzato nell'analisi dei CG

t_tester <- function(array1, array2) {
  
  difference<-sd(mapply('-', array1, array2, SIMPLIFY = T),na.rm = T)
  if ( difference != 0 || is.na(difference)) {
    #risolve l'errore Error in t.test not enough 'y' observations
    tryCatch({
      
      A <- ks.test(array1, array2)[c('p.value')]
    },error = function(error_condition){
      A<-list(p.value=NA)
    })
    
    return(c(A$p.value))
    
  }
}


Analisi <- function(matrix1){
  
  medianUP <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "UP"),-ncol(matrix1)]), 2, median,na.rm=T)
  medianMID <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "Medium"),-ncol(matrix1)]), 2, median,na.rm=T)
  medianDOWN <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "Down"),-ncol(matrix1)]), 2, median,na.rm=T)
  #inserisco alla fine della colonna di ogni gene la MEDIA DEL GRUPPO DOWN perché è sempre dispari
  df <- data.frame()
  
  df <- rbind(df, medianDOWN,medianMID,medianUP)
  
  if(dim(matrix1)[2] == 2){
    colnames(df)<-"value"
    colnames(matrix1)<-c("value","stratification")
    colnamesDfTtest<-"value"
    dimM <- 1
  }else{
    colnames(df) <- colnames(matrix1[,-ncol(matrix1)])
    dimM <- dim(matrix1[, -ncol(matrix1)])[2]
    colnamesDfTtest <- colnames(matrix1[,-ncol(matrix1)])
  }
  
  dfTtest<- setNames(data.frame(matrix(nrow = 3, ncol = dimM)), colnamesDfTtest)
  
  matrix1 <- rbind.fill(matrix1, df)
  #remove(df)
  remove(medianUP, medianMID, medianDOWN)
  
  #matrix1<- calcFC(matrix1,T)
  matrix1 <- calcBetaDifference(matrix1)
  
  for (k in 1:dimM) {
    #UPvsMID
    dfTtest[1, k]<-t_tester(matrix1[which(matrix1$stratification %in% "UP"), k], matrix1[which(matrix1$stratification %in% "Medium"), k])
    
    #UPvsDOWN
    dfTtest[2, k]<-t_tester(matrix1[which(matrix1$stratification %in% "UP"), k], matrix1[which(matrix1$stratification %in% "Down"), k])
    
    #MIDvsDOWN
    dfTtest[3, k]<-t_tester(matrix1[which(matrix1$stratification %in% "Medium"), k], matrix1[which(matrix1$stratification %in% "Down"), k])
  }
  
  matrix1 <- rbind.fill(matrix1, dfTtest)
  remove(dfTtest)
  return(matrix1)
  
}

#utilizzata per i raggruppamenti cg per isole e posizioni
Analisi2<-function(leng,index,position,column,nome_col,flag){
  
  mFinale <- data.frame(matrix())
  
  for (k in 1:leng) {
    cg <-
      as.vector(dfCGunique[which(dfCGunique$gene %in% geni[index] &
                                   dfCGunique[,column] %in% position[k]), 1])
    if (length(cg) != 0) {
      cg <- paste(cg, geni[index], sep = "_")
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
        
        m2 <- Analisi(m2)
        m2 <- as.data.frame(m2[,-2])
        
        nome_col <-
          paste(nome_col,
                paste(position[k], geni[index], sep = "_"),
                sep = ",")
        
        m2<- rbind(m2,cor(dfPancan2[c(1:473),geni[index]],m2[c(1:473),]))
        m2<- rbind(m2,cor(dfPancan2[c(1:473),geni[index]],m2[c(1:473),],method = "spearman"))
        
        m2 <- data.frame(sapply(m2, c, unlist(valExprGene[, geni[index]])), row.names = NULL)
        
        m2 <- as.data.frame(m2[-c(1:num_row),])
        mFinale <- cbind(mFinale, m2)
      }
    }
  }

  if(flag)
    assign('nome_colonne_position', nome_col, envir = .GlobalEnv)
  else
    assign('nome_colonne_island', nome_col, envir = .GlobalEnv)
  
  return(mFinale[, -1])
}

myNormShapiro <- function(x) {
  
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (is.na(n) || n < 3L || n > 5000L)
    stop("sample size must be between 3 and 5000")
  rng <- x[n] - x[1L]
  if (rng != 0) {
    if (rng < 1e-10)
      x <- x / rng
    res <- .Call(C_SWilk, x)
  }
  else{
    res <- c(NA, NA)
  }
  
  RVAL <- list(
    statistic = c(W = res[1]),
    p.value = res[2],
    method = "Shapiro-Wilk normality test",
    data.name = DNAME
  )
  class(RVAL) <- "htest"
  return(RVAL)
  
}

setRowNames <- function(df){
  row.names(df)[1] <- "medianDown"
  row.names(df)[2] <- "medianMedium"
  row.names(df)[3] <- "medianUP"
  row.names(df)[4] <- "bd_UPvsMID"
  row.names(df)[5] <- "bd_UPvsDOWN"
  row.names(df)[6] <- "bd_MIDvsDOWN"
  row.names(df)[7] <- "pvalue_UPvsMID"
  row.names(df)[8] <- "pvalue_UPvsDOWN"
  row.names(df)[9] <- "pvalue_MIDvsDOWN"
  row.names(df)[10] <- "pearson_correlation"
  row.names(df)[11] <- "spearman_correlation"
  row.names(df)[12] <- "fc_UPvsMID(gene)"
  row.names(df)[13] <- "fc_UPvsDOWN(gene)"
  row.names(df)[14] <- "fc_MIDvsDOWN(gene)"
  row.names(df)[15] <- "pvalue_UPvsMID(gene)"
  row.names(df)[16] <- "pvalue_UPvsDOWN(gene)"
  row.names(df)[17] <- "pvalue_MIDvsDOWN(gene)"
  return(df)
}