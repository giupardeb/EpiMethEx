library(readxl)
nameFolderDest <- ""
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
  )) >= 0.1 &
    (abs(
      unificati$bd_UPvsDOWN
    )) >= 0.1 & (abs(
      unificati$bd_MIDvsDOWN
    )) >= 0.1)), ]

unificati_new <-
  unificati_new[!((unificati$pvalue_UPvsMID <= 0.01) &
                    (unificati$pvalue_UPvsDOWN <= 0.01) &
                    (unificati$pvalue_MIDvsDOWN <= 0.01)
  ), ]

unificati_new <-
  unificati_new[!((unificati$fc_UPvsMID.gene. >= 2) &
                    (unificati$fc_UPvsDOWN.gene. >= 2) &
                    (unificati$fc_MIDvsDOWN.gene. >= 2)
  ), ]

unificati_new <-
  unificati_new[!((unificati$pvalue_UPvsMID.gene. <= 0.01) &
                    (unificati$pvalue_UPvsDOWN.gene. <= 0.01) &
                    (unificati$pvalue_MIDvsDOWN.gene. <= 0.01)
  ), ]

unificati_new <-
  unificati_new[!(unificati$pvalue_pearson_correlation <= 0.05), ]

write.xlsx(
  unificati_new,
  paste(nameFolderDest, "filtered_Table.xlsx", sep = "/") ,
  sheetName = "Sheet1"
)