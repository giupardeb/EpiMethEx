calcBetaDifference <- function(matrix3) {
  
  medianUP <- matrix3[nrow(matrix3), -ncol(matrix3)]
  medianMID <- matrix3[nrow(matrix3) - 1, -ncol(matrix3)]
  medianDOWN <- matrix3[nrow(matrix3) - 2, -ncol(matrix3)]
  
  if (dim(matrix3)[2] == 2) {
    bd_UPvsMID <- data.frame(value = medianUP - medianMID)
    bd_UPvsDOWN <- data.frame(value = medianUP - medianDOWN)
    bd_MIDvsDOWN <- data.frame(value = medianMID - medianDOWN)
  } else{
    bd_UPvsMID <- data.frame(medianUP - medianMID)
    bd_UPvsDOWN <- data.frame(medianUP - medianDOWN)
    bd_MIDvsDOWN <- data.frame(medianMID - medianDOWN)
    colnames(bd_UPvsMID) <- colnames(matrix3[, -ncol(matrix3)])
    colnames(bd_UPvsDOWN) <- colnames(matrix3[, -ncol(matrix3)])
    colnames(bd_MIDvsDOWN) <- colnames(matrix3[, -ncol(matrix3)])
  }
  
  matrix3 <- rbind.fill(matrix3, bd_UPvsMID, bd_UPvsDOWN, bd_MIDvsDOWN)
  
  return(matrix3)
}

calcFC <- function(matrix2, dataIsLinear) {
  
  df <- data.frame()
  
  meanUp <- matrix2[nrow(matrix2), -ncol(matrix2)]
  meanMid <- matrix2[nrow(matrix2) - 1, -ncol(matrix2)]
  meanDown <- matrix2[nrow(matrix2) - 2, -ncol(matrix2)]
  
  if (dataIsLinear) {
    fc_UPvsMID <- setFC(meanUp, meanMid)
    fc_UPvsDOWN <-  setFC(meanUp, meanDown)
    fc_MIDvsDOWN <- setFC(meanMid, meanDown)
  } else{
    fc_UPvsMID <- setFClog(meanUp, meanMid)
    fc_UPvsDOWN <- setFClog(meanUp, meanDown)
    fc_MIDvsDOWN <-  setFClog(meanMid, meanDown)
  }
  
  df <- rbind(df, fc_UPvsMID, fc_UPvsDOWN, fc_MIDvsDOWN)
  colnames(df) <- colnames(matrix2[,-ncol(matrix2)])
  
  matrix2 <- rbind.fill(matrix2, df)
  remove(df)
  
  return(matrix2)
}


checkSign <- function(a, b) {
  
  return (sign(a) == sign(b))
}

setFClog <- function(meanFirstGroup, meanSecondGroup) {
  
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

#Used in genes analysis
ttester <- function(array1, array2) {
  
  if (sd(mapply('-', array1, array2, SIMPLIFY = T), na.rm = T) != 0) {
    A <- t.test(array1, array2,
                var.equal = F)[c('statistic', 'p.value')]
    return(c(A$statistic, A$p.value))
  } else{
    return(c(NA, NA))
  }
}

t_tester <- function(array1, array2) {
  
  difference <- sd(mapply('-', array1, array2, SIMPLIFY = T), na.rm = T)
  
  if (difference != 0 || is.na(difference)) {
    tryCatch({
      A <- ks.test(array1, array2)[c('p.value')]
    }, error = function(error_condition) {
      A <- list(p.value = NA)
    })
    
    return(c(A$p.value))
  }
}

Analisi <- function(matrix1) {
  
  medianUP <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "UP"),-ncol(matrix1)]), 2, median, na.rm = T)
  medianMID <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "Medium"),-ncol(matrix1)]), 2, median, na.rm = T)
  medianDOWN <- apply(as.data.frame(matrix1[which(matrix1$stratification %in% "Down"),-ncol(matrix1)]), 2, median, na.rm = T)

  if (dim(matrix1)[2] == 2) {
    colnames(matrix1) <- c("value", "stratification")
    colnamesDfTtest <- "value"
    dimM <- 1
  } else{
    dimM <- dim(matrix1[, -ncol(matrix1)])[2]
    colnamesDfTtest <- colnames(matrix1[,-ncol(matrix1)])
  }
  
  dfTtest <- setNames(data.frame(matrix(nrow = 3, ncol = dimM)), colnamesDfTtest)
  
  matrix1 <- rbind(matrix1, medianDOWN, medianMID, medianUP)
  
  remove(medianUP, medianMID, medianDOWN)
  
  matrix1 <- calcBetaDifference(matrix1)
  
  for (k in 1:dimM) {
    #UPvsMID
    dfTtest[1, k] <-
      t_tester(matrix1[which(matrix1$stratification %in% "UP"), k], matrix1[which(matrix1$stratification %in% "Medium"), k])
    #UPvsDOWN
    dfTtest[2, k] <-
      t_tester(matrix1[which(matrix1$stratification %in% "UP"), k], matrix1[which(matrix1$stratification %in% "Down"), k])
    #MIDvsDOWN
    dfTtest[3, k] <-
      t_tester(matrix1[which(matrix1$stratification %in% "Medium"), k], matrix1[which(matrix1$stratification %in% "Down"), k])
  }
  
  matrix1 <- rbind.fill(matrix1, dfTtest)
  remove(dfTtest)
  return(matrix1)
}

#Used for islands and positions cg groupings
Analisi2 <- function(leng, index, position, column) {
  
  mFinale <- data.frame(matrix())
  for (k in 1:leng) {
    cg <- as.vector(dfCGunique[which(dfCGunique$gene %in% genes[index] & dfCGunique[, column] %in% position[k]), 1])
    
    if (length(cg) != 0) {
      cg <- paste(cg, genes[index], sep = "_")
      m2 <- as.data.frame(m[, cg])
      m2 <- as.data.frame(m2[, colSums(is.na(m2)) != nrow(m2)]) #remove all coloumns that have all values NA
      
      if (dim(m2)[2] > 1) {
        m2 <- stack(m2)
        m2 <- as.matrix(m2[,-2])
        num_row <- nrow(m2)
        m2 <- cbind(m2, stratification)
        colnames(m2) <- c("value", "stratification")
        
        m2 <- Analisi(m2)
        m2 <- as.data.frame(m2[,-2])
        
        c <- as.data.frame(rep(dfPancan2[c(1:473), genes[i]], length(cg)))
        mTmp <- as.data.frame(m2[1:dim(c)[1], ])
        
        a <- corr.test(c, mTmp, adjust = "none")
        m2 <- rbind(m2, as.numeric(a$r), as.numeric(a$p)) #add pearson correlation and p-value
        
        m2 <- data.frame(sapply(m2, c, unlist(valExprGene[, genes[index]])), row.names = NULL)
        
        m2 <- as.data.frame(m2[-c(1:num_row),])
        colnames(m2) <- paste(position[k], genes[index], sep = "_")
        mFinale <- cbind(mFinale, m2)
      }
    }
  }
  return(subset(mFinale, select = -c(1)))
}

setRowNames <- function(df) {
  
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
  row.names(df)[11] <- "pvalue_pearson_correlation"
  row.names(df)[12] <- "fc_UPvsMID(gene)"
  row.names(df)[13] <- "fc_UPvsDOWN(gene)"
  row.names(df)[14] <- "fc_MIDvsDOWN(gene)"
  row.names(df)[15] <- "pvalue_UPvsMID(gene)"
  row.names(df)[16] <- "pvalue_UPvsDOWN(gene)"
  row.names(df)[17] <- "pvalue_MIDvsDOWN(gene)"
  return(df)
}
