ordinamento <- function(df, nameDF) {
  A <- as.data.frame(t(df))
  colnames(A) <- c("FC", "pValue")
  A <- cbind(A, V22 = rownames(A))
  
  if (nrow(A) > n) {
    A <- A[order(A$FC, decreasing = T),]
    A <- rbind(head(A, n / 2), tail(A, n / 2))
  }
  
  C <- merge(A, righeCheTiServono1, sort = F)
  
  colnames(C) <-
    c("Gene", "FC", "pValue", "CG", "Identificativo", "Posizione")
  
  #raggruppa FC GENE PVALUE, i cg sono all'interno di posizione
  Z <- C %>%
    group_by(Gene = C$Gene,
             FC = C$FC,
             pValue = C$pValue) %>%
    summarise(Posizione = paste(
      CG,
      Identificativo,
      Posizione,
      sep = ";;",
      collapse = ";;"
    ))
  
  df <- Z[order(Z$FC, decreasing = T),]
  
  assign(nameDF, df, envir = .GlobalEnv)
}