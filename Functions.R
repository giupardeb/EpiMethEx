calcFC <- function(matrix2,dataIsLinear) {
  df <- data.frame()
  
  meanUp <- matrix2[476, -ncol(matrix2)]
  meanMid <- matrix2[475, -ncol(matrix2)]
  meanDown <- matrix2[474, -ncol(matrix2)]
  
  if(dataIsLinear){
    fc_UPvsMID <- setFC(meanUp, meanMid)
    fc_UPvsDOWN <-  setFC(meanUp, meanDown)
    fc_MIDvsDOWN <- setFC(meanMid, meanDown)
  }else{
    fc_UPvsMID <- setFClog(meanUp,meanMid) 
    fc_UPvsDOWN <- setFClog(meanUp,meanDown)
    fc_MIDvsDOWN <-  setFClog(meanMid,meanDown)
  }
  
  df <- rbind(df, fc_UPvsMID)
  df <- rbind(df, fc_UPvsDOWN)
  df <- rbind(df, fc_MIDvsDOWN)
  
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
ttester <- function(array1, array2, start, end) {
  
  if (sd(mapply('-', array1, array2, SIMPLIFY = T)) != 0) {
    A <- t.test(array1, array2,
                var.equal = F)[c('statistic', 'p.value')]
    dfTtest[start:end, k] <- c(A$statistic, A$p.value)
    
    assign('dfTtest', dfTtest, envir = .GlobalEnv)
  }
  
}

#Utilizzato nell'analisi dei CG

t_tester <- function(array1, array2) {
  
  if (sd(mapply('-', array1, array2, SIMPLIFY = T),na.rm = T) != 0) {
    
    A <- t.test(array1, array2,
                var.equal = F)[c('p.value')]
    
    return(c(A$p.value))
    
  }
}


Analisi <- function(matrix1){
  
  meanUP <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "UP"),-ncol(matrix1)]), 2, mean,na.rm=T)
  meanMID <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "Medium"),-ncol(matrix1)]), 2, mean,na.rm=T)
  meanDOWN <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "Down"),-ncol(matrix1)]), 2, mean,na.rm=T)
  #inserisco alla fine della colonna di ogni gene la MEDIA DEL GRUPPO DOWN perché è sempre dispari
  df <- data.frame()
  
  df <- rbind(df, meanDOWN)
  df <- rbind(df, meanMID)
  df <- rbind(df, meanUP)
  
  rownames(df) <- c("meanDOWN", "meanMID", "meanUP")
  
  if(dim(matrix1)[2] == 2){
    colnames(df)<-"value"
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
  remove(meanUP, meanMID, meanDOWN)
  
  matrix1<- calcFC(matrix1,T)
  
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
