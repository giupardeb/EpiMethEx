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
        
        #eliminare le righe dei valori dei cg
        m2 <- as.data.frame(m2[-c(1:num_row),])
        m2 <- data.frame(sapply(m2, c, unlist(valExprGene[, geni[index]])), row.names = NULL)
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