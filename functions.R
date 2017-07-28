calcFCGeneLinear <- function() {
  df <- data.frame()
  
  meanUp <- dfPancan2[476, ]
  meanMid <- dfPancan2[475, ]
  meanDown <- dfPancan2[474, ]
  
  fc_UPvsMID <- setFC(meanUp, meanMid)
  fc_UPvsDOWN <-  setFC(meanUp, meanDown)
  fc_MIDvsDOWN <- setFC(meanMid, meanDown)
  
  df <- rbind(df, fc_UPvsMID)
  df <- rbind(df, fc_UPvsDOWN)
  df <- rbind(df, fc_MIDvsDOWN)
  rownames(df) <- c("fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN")
  colnames(df) <- colnames(dfPancan2)
  dfPancan2 <- rbind(dfPancan2, df)
  remove(df)
  assign('dfPancan2', dfPancan2, envir = .GlobalEnv)
}

calcFCLinear <- function(a, b) {
  #df <- data.frame()
  
  #meanUp <- dfPancan2[476,]
  #meanMid <- dfPancan2[475,]
  #meanDown <- dfPancan2[474,]
  
  fc <- setFC(a, b)
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
  
  df <- data.frame()
  
  meanUp <- dfPancan2[476, ]
  meanMid <- dfPancan2[475, ]
  meanDown <- dfPancan2[474, ]
  
  fc_UPvsMID <- setFClog(meanUp,meanMid ) 
  fc_UPvsDOWN <- setFClog(meanUp,meanDown)
  fc_MIDvsDOWN <-  setFClog(meanMid,meanDown)
  
  df <- rbind(df, fc_UPvsMID)
  df <- rbind(df, fc_UPvsDOWN)
  df <- rbind(df, fc_MIDvsDOWN)
  rownames(df) <- c("fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN")
  names(df)<-names(dfPancan2)
  dfPancan2 <- rbind(dfPancan2, df)
  remove(df)
  assign('dfPancan2', dfPancan2, envir = .GlobalEnv)
  
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

ttester <- function(array1, array2, start, end) {
  
  if (sd(mapply('-', array1, array2, SIMPLIFY = T)) != 0) {
    A <- t.test(array1, array2,
                var.equal = F)[c('statistic', 'p.value')]
    dfTtest[start:end, k] <- c(A$statistic, A$p.value)
    
    assign('dfTtest', dfTtest, envir = .GlobalEnv)
  }
  
}