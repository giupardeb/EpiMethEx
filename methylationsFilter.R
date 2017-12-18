library(readxl)
nameFolderDest <- ""
Pvalue <- 0.01
FoldChange <- 2
PvaluePearson <- 0.05
BetaDifference <- 0.1

unificati <-
  read_excel(file.choose())

unificati_new <-
  unificati[!(((unificati$medianUP > unificati$medianMedium) &
                 (unificati$medianMedium > unificati$medianDown)
  ) |
    ((unificati$medianUP < unificati$medianMedium) &
       (unificati$medianMedium < unificati$medianDown)
    )), ]

unificati_new <-
  unificati_new[!(((abs(
    unificati$bd_UPvsMID
  )) >= BetaDifference &
    (abs(
      unificati$bd_UPvsDOWN
    )) >= BetaDifference & (abs(
      unificati$bd_MIDvsDOWN
    )) >= BetaDifference)), ]

unificati_new <-
  unificati_new[!((unificati$pvalue_UPvsMID <= Pvalue) &
                    (unificati$pvalue_UPvsDOWN <= Pvalue) &
                    (unificati$pvalue_MIDvsDOWN <= Pvalue)
  ), ]

unificati_new <-
  unificati_new[!((unificati$fc_UPvsMID.gene. >= FoldChange) &
                    (unificati$fc_UPvsDOWN.gene. >= FoldChange) &
                    (unificati$fc_MIDvsDOWN.gene. >= FoldChange)
  ), ]

unificati_new <-
  unificati_new[!((unificati$pvalue_UPvsMID.gene. <= Pvalue) &
                    (unificati$pvalue_UPvsDOWN.gene. <= Pvalue) &
                    (unificati$pvalue_MIDvsDOWN.gene. <= Pvalue)
  ), ]

unificati_new <-
  unificati_new[!(unificati$pvalue_pearson_correlation <= PvaluePearson), ]

write.xlsx(
  unificati_new,
  paste(nameFolderDest, "filtered_Table.xlsx", sep = "/") ,
  sheetName = "Sheet1"
)
