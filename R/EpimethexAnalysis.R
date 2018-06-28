#' Create a Epimethex function
#'
#' @param path_dfPancan String, is the path of expression data
#' @param path_dfGpl String, is the path of Annotations data
#' @param path_dfMethylation String, is the path of Methylation data
#' @param minRangeGene Numeric
#' @param maxRangeGene Numeric
#' @param numCores Numeric, is the number of cores that you can use
#' @param dataLinear Logic.
#' @param outputPath String
#'
#' @import doParallel
#' @import foreach
#' @import miscTools
#' @import stats
#' @import utils
#' @import magrittr
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @source Functions
#'
#' @export
epimethex.analysis <- function(path_dfPancan,
                               path_dfGpl,
                               path_dfMethylation,
                               minRangeGene,
                               maxRangeGene,
                               numCores,
                               dataLinear,
                               outputPath) {
  dfPancan <-
    read.table(path_dfPancan, header = T,
               sep = "\t")[minRangeGene:maxRangeGene, ]

  # name of folder where files xml will saved
  nameFD <- paste(as.character(minRangeGene), as.character(maxRangeGene), sep = "-")
  dir.create(file.path(outputPath, nameFD))

  # the first 37 rows of probset dataset, are comments.
  # The following istruction read only columns of interest
  dfGpl <-
    data.table::fread(
      path_dfGpl,
      sep = "\t",
      header = T ,
      skip = 37,
      na.strings = c("", "NA")
    )[-1, c(
      'ID',
      'Relation_to_UCSC_CpG_Island',
      'UCSC_CpG_Islands_Name',
      'UCSC_RefGene_Accession',
      'Chromosome_36',
      'Coordinate_36',
      'UCSC_RefGene_Name',
      'UCSC_RefGene_Group'
    )]

  # create a new column `island` with the two columns collapsed together
  dfGpl$island <-
    apply(dfGpl[, c('Relation_to_UCSC_CpG_Island', 'UCSC_CpG_Islands_Name')] , 1 , paste , collapse = "_")

  # remove the unnecessary columns

  dfGpl <-
    dfGpl[,-which(names(dfGpl) %in% c('Relation_to_UCSC_CpG_Island', 'UCSC_CpG_Islands_Name'))]

  dfMethylation <-
    read.table(path_dfMethylation,
               sep = "\t",
               header = T)

  # transpose all but the first column (name)
  dfPancan2 <- as.data.frame(t(dfPancan[, -1]))
  colnames(dfPancan2) <- dfPancan$sample

  #remove all genes that has all values equal a zero
  dfPancan2 <- dfPancan2[, which(!apply(dfPancan2 == 0, 2, all))]

  #indexTCGA contains gene by gene the patient decreasing ordered
  indexTCGA <-
    data.frame(sapply(seq_along(dfPancan2), function(x) {
      row.names(dfPancan2[order(dfPancan2[x], decreasing = T), x, drop = F])
    }))
  colnames(indexTCGA) <- colnames(dfPancan2)

  dfPancan2 <- apply(dfPancan2, 2, sort, decreasing = T)

  #calculate quantili of all genes into dfPancan2
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    quantili <-
      data.frame(matrixStats::colQuantiles(dfPancan2, probs = c(0.33, 0.66, 0.99)))
  }

  data.table::setnames(quantili, c("perc33", "perc66", "perc99"))
  quantili <- as.data.frame(t(quantili))

  dfPancan2 <- as.data.frame(dfPancan2)

  # [START] these istruction divide into equal parts the dfPancan2
  index <-
    which.max(apply(dfPancan2, 2, function(x)
      length(unique(x))))[[1]]

  #if (requireNamespace("RcmdrMisc", quietly = TRUE)) {
  dfPancan2$variable <-
    with(
      dfPancan2,
      RcmdrMisc::bin.var(
        dfPancan2[, index],
        bins = 3,
        method = 'proportions',
        labels = c('D', 'M', 'UP')
      )
    )
  #}

  # [END]
  lengthDfpancan <- nrow(dfPancan2)
  # [START] calculate the means of all values of genes that are UP, medium and down
  meanUP <-
    apply(dfPancan2[which(dfPancan2$variable %in% "UP"), -ncol(dfPancan2)], 2, mean)
  meanMID <-
    apply(dfPancan2[which(dfPancan2$variable %in% "M"), -ncol(dfPancan2)], 2, mean)
  meanDOWN <-
    apply(dfPancan2[which(dfPancan2$variable %in% "D"), -ncol(dfPancan2)], 2, mean)
  # [END]
  dfPancan2 <- rbind(dfPancan2, meanDOWN, meanMID, meanUP)

  remove(meanUP, meanMID, meanDOWN)

  dfPancan2 <- plyr::rbind.fill(dfPancan2, quantili)

  remove(quantili)

  ROW_MEAN_DOWN <- nrow(dfPancan2) - 5
  ROW_MEAN_MEDIUM <- nrow(dfPancan2) - 4
  ROW_MEAN_UP <- nrow(dfPancan2) - 3
  ROW_PERC_33 <- nrow(dfPancan2) - 2
  ROW_PERC_66 <- nrow(dfPancan2) - 1
  ROW_PERC_99 <- nrow(dfPancan2)
  gc()

  #Calculate Fold change of the combinations (up vs mid, up vs down, etc..)
  dfPancan2 <- calcFC(dfPancan2, dataLinear)

  dimDFpancan <- dim(dfPancan2[,-ncol(dfPancan2)])[2]

  # Create cluster with desired number of cores
  cl <- parallel::makeCluster(numCores)

  # Register cluster
  registerDoParallel(cl)

  #[START] calculate t-test for all genes
  resultUPvsMID <-
    foreach(indexTmp = 1:dimDFpancan, .combine = rbind) %dopar% {
      #UPvsMID
      calculateTtest(dfPancan2[which(dfPancan2$variable %in% "UP"), indexTmp],
                     dfPancan2[which(dfPancan2$variable %in% "M"), indexTmp], T)
    }
  resultUPvsDOWN <-
    foreach(indexTmp = 1:dimDFpancan, .combine = rbind) %dopar% {
      #UPvsDOWN
      calculateTtest(dfPancan2[which(dfPancan2$variable %in% "UP"), indexTmp],
                     dfPancan2[which(dfPancan2$variable %in% "D"), indexTmp], T)

    }

  resultMIDvsDOWN <-
    foreach(indexTmp = 1:dimDFpancan, .combine = rbind) %dopar% {
      #MIDvsDOWN
      calculateTtest(dfPancan2[which(dfPancan2$variable %in% "M"), indexTmp],
                     dfPancan2[which(dfPancan2$variable %in% "D"), indexTmp], T)
    }
  parallel::stopCluster(cl)
  #[END]

  resultUPvsMID <- as.data.frame(t(resultUPvsMID))
  resultUPvsDOWN <- as.data.frame(t(resultUPvsDOWN))
  resultMIDvsDOWN <- as.data.frame(t(resultMIDvsDOWN))
  colnames(resultUPvsMID) <-
    colnames(dfPancan2[,-ncol(dfPancan2)])
  colnames(resultUPvsDOWN) <-
    colnames(dfPancan2[,-ncol(dfPancan2)])
  colnames(resultMIDvsDOWN) <-
    colnames(dfPancan2[,-ncol(dfPancan2)])

  dfPancan2 <-
    plyr::rbind.fill(dfPancan2, resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)

  remove(resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
  gc()

  ROW_FC_UPvsMID <- nrow(dfPancan2) - 8
  ROW_FC_UPvsDOWN <- nrow(dfPancan2) - 7
  ROW_FC_MIDvsDOWN <- nrow(dfPancan2) - 6
  ROW_TTEST_UPvsMID <- nrow(dfPancan2) - 5
  ROW_PVALUE_UPvsMID <- nrow(dfPancan2) - 4
  ROW_TTEST_UPvsDOWN <- nrow(dfPancan2) - 3
  ROW_PVALUE_UPvsDOWN <- nrow(dfPancan2) - 2
  ROW_TTEST_MIDvsDOWN <- nrow(dfPancan2) - 1
  ROW_PVALUE_MIDvsDOWN <- nrow(dfPancan2)

  rowNames <-
    c(
      "meanDown",
      "meanMedium",
      "meanUP",
      "perc33",
      "perc66",
      "perc99",
      "fc_UPvsMID",
      "fc_UPvsDOWN",
      "fc_MIDvsDOWN",
      "ttest_UPvsMID",
      "pvalue_UPvsMID",
      "ttest_UPvsDOWN",
      "pvalue_UPvsDOWN",
      "ttest_MIDvsDOWN",
      "pvalue_MIDvsDOWN"
    )
  j <- 1
  for (i in ROW_MEAN_DOWN:ROW_PVALUE_MIDvsDOWN) {
    row.names(dfPancan2)[i] <- rowNames[j]
    j <- j + 1
  }

  # [START] the dfAnnotations contains all cg with the corrispondent genes,
  # islands, position of body and position into chromosoma
  s <- strsplit(dfGpl$UCSC_RefGene_Name, split = ";")
  s1 <- strsplit(dfGpl$UCSC_RefGene_Accession, split = ";")
  s2 <- strsplit(dfGpl$UCSC_RefGene_Group, split = ";")

  dfAnnotations <-
    data.frame(
      cg = rep(dfGpl$ID, sapply(s, length)),
      #V1
      Chromosome_36 = rep(dfGpl$Chromosome_36, sapply(s, length)),
      Coordinate_36 = rep(dfGpl$Coordinate_36, sapply(s, length)),
      gene = unlist(s),
      #V22
      UCSC_RefGene_Accession = unlist(s1),
      position = unlist(s2),
      #V24
      island = rep(dfGpl$island, sapply(s, length))
    )
  #[END]

  remove(s, s1, s2, dfGpl)
  gc()

  #remove genes and cg aren't into dfPancan2 from dfAnnotations
  dfAnnotations <-
    subset(dfAnnotations,
           as.character(dfAnnotations$gene) %in% names(dfPancan2))

  #remove genes and cg aren't into dfPancan2 from dfMethylation
  dfMethylation <-
    subset(dfMethylation,
           as.character(dfMethylation$sample) %in% dfAnnotations$cg)

  maxOccurence <-
    max(as.data.frame(table(unlist(
      dfAnnotations$gene
    )))$Freq)

  #reorder genes related CG
  dfAnnotations <-
    dfAnnotations %>% plyr::arrange(dfAnnotations$gene, dfAnnotations$Coordinate_36)

  #dfCGunique contain the unique couple (CG-position)

  dfCGunique <-
    dfAnnotations[!duplicated(dfAnnotations[, c(1, 6)]), ]

  dfMethylation$sample <- as.factor(dfMethylation$sample)
  dfMethylation <- as.data.frame(dfMethylation)
  rownames(dfMethylation) <- dfMethylation$sample
  data.table::setkey(data.table::as.data.table(dfMethylation), sample)
  gc()

  tmp <- tapply(lapply(1:nrow(dfCGunique),
                       function (i)
                         (dfMethylation[as.vector(dfCGunique[i, 1]),
                                        as.vector(indexTCGA[, as.vector(dfCGunique[i, "gene"])])])),
                factor(dfCGunique[, "gene"]), function (x) {
                  unname(unlist(x))
                })

  max.rows <- max(sapply(tmp, length))

  # DFCGorder is a dataframe that contains all CG values
  # ordered by indexTCGA dataframe
  DFCGorder <-
    do.call(cbind, lapply(tmp, function(x) {
      length(x) <- max.rows
      return (x)
    }))

  valExprGene <-
    dfPancan2[c(
      ROW_FC_UPvsMID:ROW_FC_MIDvsDOWN,
      ROW_PVALUE_UPvsMID,
      ROW_PVALUE_UPvsDOWN,
      ROW_PVALUE_MIDvsDOWN
    ),-ncol(dfPancan2)]
  positions <- as.vector(unique(dfCGunique$position))
  genes <- colnames(DFCGorder)
  islands <- as.vector(unique(dfCGunique$island))
  islands <-
    Filter(function(x)
      ! any(grepl("NA_NA", x)), islands) #remove islands that have values "NA_NA"

  stratification <-
    as.data.frame(rep(c("UP", "Medium", "Down"), each = 158))
  colnames(stratification) <- "stratification"

  gc()
  lengthGens <- length(genes)
  lengthPositions <- length(positions)
  lengthIslands <- length(islands)

  # Create cluster with desired number of cores
  cl <- parallel::makeCluster(numCores)

  # Register cluster
  registerDoParallel(cl)

  parallel::clusterCall(cl,
                        requireNamespace,
                        package = "miscTools",
                        quietly = TRUE)
  parallel::clusterCall(cl,
                        requireNamespace,
                        package = "psych",
                        quietly = TRUE)
  parallel::clusterCall(cl,
                        requireNamespace,
                        package = "plyr",
                        quietly = TRUE)

  # [START] in questo ciclo for verrà creato un dataframe,
  # in cui all'interno di esso ci saranno tutte le medie,le beta difference
  # e i pvalue delle varie stratificazioni(UP,MED,DOWN), d
  # i tutti i CG in tutte le posizioni di tutti i geni.

  mFinaleCGglobali <-
    foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
      flag <- F

      # [START] creo la matrice tempMatrix che conterrà i valori
      # dei cg corrispondenti all'iesimo gene
      tempMatrix <-
        data.frame(matrix(
          DFCGorder[, i],
          nrow = lengthDfpancan,
          ncol = length(DFCGorder[, i]) / lengthDfpancan
        ))
      colnames(tempMatrix) <-
        as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
      #[END]

      # [START] elimino in tempMatrix le colonne che non sono
      # associate a nessun cg dell'iesimo gene
      keep.cols <- names(tempMatrix) %in% NA
      tempMatrix <- tempMatrix [!keep.cols]
      #[END]

      colnames(tempMatrix) <-
        paste(names(tempMatrix), genes[i], sep = "_")

      #conservo su columnNA gli indici delle colonne che hanno tutti i valori NA
      columnNA <- which(sapply(tempMatrix, function(x)
        all(is.na(x))))

      tempMatrix <-
        tempMatrix[, colSums(is.na(tempMatrix)) != nrow(tempMatrix)]

      if (length(tempMatrix) != 0) {
        tempMatrix <- cbind(tempMatrix, stratification)
        tempMatrix <- Analysis(tempMatrix)
        tempMatrix$stratification <- NULL
      } else{
        flag <- T
      }

      # inserisco in tempMatrix le colonne che precedentemente
      # ho eliminato nella stessa posizione

      if (length(columnNA) != 0) {
        for (j in 1:length(columnNA)) {
          tempMatrix <-
            insertCol(as.matrix(tempMatrix), columnNA[[j]], v = NA)
        }
      }

      if (flag) {
        df <- as.data.frame(matrix(NA, nrow = 9, ncol = dim(tempMatrix)[2]))
        tempMatrix <- rbind(tempMatrix, df)
      }

      # [START] calcolo la correlazione tra i dati di espressione
      # dell'i-esimo gene con i dati di metilazione dei cg associati
      DataexpressionGeneTmp <-
        as.data.frame(dfPancan2[c(1:lengthDfpancan), genes[i]])
      DataCG_Tmp <-
        as.data.frame(tempMatrix[c(1:lengthDfpancan),])
      # calcolo il test di correlazione tra i dati di espressione
      # del gene e con i relativi CG
      resultCorrTest <-
        psych::corr.test(DataexpressionGeneTmp, DataCG_Tmp , adjust = "none")

      tempMatrix <-
        rbind(tempMatrix,
              as.numeric(resultCorrTest$r),
              as.numeric(resultCorrTest$p))

      # [START] elimino da tempMatrix i valori di espressione dei cg, tenendo
      # soltando le righe: "medianDown", "medianMedium", "medianUP",
      # "bd_UPvsMID", "bd_UPvsDOWN", "bd_MIDvsDOWN", "pvalue_UPvsMID",
      # "pvalue_UPvsDOWN", "pvalue_MIDvsDOWN", "pearson_correlation",
      #  "pvalue_pearson_correlation"

      if (dim(tempMatrix)[1] > lengthDfpancan) {
        tempMatrix <- as.data.frame(tempMatrix[-c(1:lengthDfpancan),])
      } else{
        tempMatrix <- as.data.frame(tempMatrix[-c(1:464),])
      }
      #[END]

      # insert into temp matrix "tempMatrix" the following values of the
      # gene[i]: "fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN",
      # "pvalue_UPvsMID", "pvalue_UPvsDOWN", "pvalue_MIDvsDOWN"
      tempMatrix <-
        data.frame(sapply(tempMatrix, c, unlist(valExprGene[, genes[i]])), row.names = NULL)

      # [START] cambio il nome delle colonne associacndo il nome del cg
      # con il corrispettivo posizione e gene. ex. cg27394127_TSS1500_BFAR
      colnames(tempMatrix) <-
        paste(as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1]),
              as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 6]),
              genes[i],
              sep = "_")
      #[END]
      remove(DataexpressionGeneTmp, DataCG_Tmp, resultCorrTest)

      #print the matrix1 into mFinaleCGglobali
      tempMatrix
    }
  #[END]

  mFinaleCGglobali <- setRowNames(mFinaleCGglobali)

  # [START] estraggo i nomi dei cg e successivamente vado a
  # recuperare tutte le isole associate a quei cg.
  tmp <- strsplit(colnames(mFinaleCGglobali), split = "_")
  tmp <- unlist(lapply(tmp, '[[', 1))

  mFinaleCGglobali <-
    rbind(mFinaleCGglobali, as.character(dfCGunique[which(dfCGunique$cg %in% as.array(tmp)), 7]))
  row.names(mFinaleCGglobali)[nrow(mFinaleCGglobali)] <- "island"
  #[END]

  mFinaleCGglobali <- t(mFinaleCGglobali)

  xlsx::write.xlsx(
    mFinaleCGglobali,
    file.path(outputPath, nameFD, "CG_data_individually.xlsx"),
    sheetName = "Sheet1"
  )
  remove(mFinaleCGglobali, tmp)
  gc()

  mFinaleCGunificati <-
    foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
      # [START] creo la matrice tempMatrix che conterrà i valori dei cg
      # corrispondenti all'iesimo gene
      tempMatrix <-
        data.frame(matrix(
          DFCGorder[, i],
          nrow = lengthDfpancan,
          ncol = length(DFCGorder[, i]) / lengthDfpancan
        ))
      colnames(tempMatrix) <-
        as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
      #[END]

      # [START] elimino in tempMatrix le colonne che non sono
      # associate a nessun cg dell'iesimo gene
      keep.cols <- names(tempMatrix) %in% NA
      tempMatrix <- tempMatrix [!keep.cols]
      #[END]

      colnames(tempMatrix) <-
        paste(names(tempMatrix), genes[i], sep = "_")

      tempMatrix <-
        as.data.frame(tempMatrix[, colSums(is.na(tempMatrix)) != nrow(tempMatrix)])
      num_CG <- length(tempMatrix)

      if (dim(tempMatrix)[2] != 0) {
        if (dim(tempMatrix)[2] != 1) {
          tempMatrix <- stack(tempMatrix)
          tempMatrix <- as.matrix(tempMatrix[, -2])
        }
        num_row_m4 <- nrow(tempMatrix)
        tempMatrix <- cbind(tempMatrix, stratification)
        colnames(tempMatrix) <- c("value", "stratification")
        tempMatrix <- Analysis(tempMatrix)
        tempMatrix <- as.data.frame(tempMatrix[, -2])

        names(tempMatrix) <- "value"

        # [START] calcolo la correlazione tra i dati di espressione
        #  dell'i-esimo gene con i dati di metilazione dei cg associati
        dfTmp <-
          as.data.frame(rep(dfPancan2[c(1:lengthDfpancan), genes[i]], num_CG))
        m4Tmp <- as.data.frame(tempMatrix[c(1:num_row_m4),])
        resultCorrTest <-
          psych::corr.test(dfTmp, m4Tmp, adjust = "none")


        tempMatrix <-
          rbind(tempMatrix,
                as.numeric(resultCorrTest$r),
                as.numeric(resultCorrTest$p))
        #[END]

        tempMatrix <-
          as.data.frame(tempMatrix[-c(1:num_row_m4), ])
        tempMatrix <-
          data.frame(sapply(tempMatrix, c, unlist(valExprGene[, genes[i]])), row.names = NULL)
        colnames(tempMatrix) <- paste("CG", genes[i], sep = "_")
        remove(dfTmp, m4Tmp, resultCorrTest)

        tempMatrix

      } else{
        tempMatrix <- as.data.frame(matrix(NA, nrow = 17, ncol = 1))
        colnames(tempMatrix) <- paste("CG", genes[i], sep = "_")
        tempMatrix
      }
    }

  mFinaleCGunificati <- setRowNames(mFinaleCGunificati)

  mFinaleCGunificati <- t(mFinaleCGunificati)

  xlsx::write.xlsx(
    mFinaleCGunificati,
    file.path(outputPath, nameFD, "CG_of_a_genes.xlsx"),
    sheetName = "Sheet1"
  )

  ###CG group by gene position
  mFinaleCGposition <-
    foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
      # [START] creo la matrice tempMatrix che conterrà i valori
      # dei cg corrispondenti all'iesimo gene
      tempMatrix <-
        data.frame(matrix(
          DFCGorder[, i],
          nrow = lengthDfpancan,
          ncol = length(DFCGorder[, i]) / lengthDfpancan
        ))
      colnames(tempMatrix) <-
        as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
      #[END]

      # [START] elimino in tempMatrix le colonne che non sono
      # associate a nessun cg dell'iesimo gene
      keep.cols <- names(tempMatrix) %in% NA
      tempMatrix <- tempMatrix [!keep.cols]
      #[END]

      colnames(tempMatrix) <-
        paste(names(tempMatrix), genes[i], sep = "_")

      AnalysisIslands_PositionsCG(
        lengthPositions,
        i,
        positions,
        "position",
        dfCGunique,
        genes,
        tempMatrix,
        stratification,
        dfPancan2,
        valExprGene
      )

    }

  mFinaleCGposition <- setRowNames(mFinaleCGposition)

  mFinaleCGposition <- t(mFinaleCGposition)

  xlsx::write.xlsx(
    mFinaleCGposition,
    file.path(outputPath, nameFD, "CG_by_position.xlsx"),
    sheetName = "Sheet1"
  )

  #CG group by island CpG
  mFinaleCGisland <-
    foreach(i = 1:lengthGens, .combine = cbind) %dopar% {
      # [START] creo la matrice tempMatrix che conterrà i
      # valori dei cg corrispondenti all'iesimo gene
      tempMatrix <-
        data.frame(matrix(
          DFCGorder[, i],
          nrow = lengthDfpancan,
          ncol = length(DFCGorder[, i]) / lengthDfpancan
        ))
      colnames(tempMatrix) <-
        as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), 1])
      #[END]

      # [START] elimino in tempMatrix le colonne che non sono
      # associate a nessun cg dell'iesimo gene
      keep.cols <- names(tempMatrix) %in% NA
      tempMatrix <- tempMatrix [!keep.cols]
      #[END]

      colnames(tempMatrix) <-
        paste(names(tempMatrix), genes[i], sep = "_")

      AnalysisIslands_PositionsCG(
        lengthIslands,
        i,
        islands,
        "island",
        dfCGunique,
        genes,
        tempMatrix,
        stratification,
        dfPancan2,
        valExprGene
      )

    }

  mFinaleCGisland <- setRowNames(mFinaleCGisland)

  mFinaleCGisland <- t(mFinaleCGisland)

  xlsx::write.xlsx(mFinaleCGisland,
                   file.path(outputPath, nameFD, "CG_by_islands.xlsx"),
                   sheetName = "Sheet1")

}
