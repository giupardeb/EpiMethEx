#' Create a Epimethex function
#'
#' @param dfPancan data.frame, expression data
#' @param dfGpl data.frame, Annotations data
#' @param dfMethylation data.frame, Methylation data
#' @param minRangeGene Numeric
#' @param maxRangeGene Numeric
#' @param numCores Numeric, is the number of cores that you can use
#' @param dataLinear Logic.
#'
#' @import doParallel
#' @import foreach
#' @import miscTools
#' @import stats
#' @import utils
#' @import magrittr
#' @import parallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @examples
#'
#' Annotations <- data.frame(
#'     ID = c("cg1","cg2", "cg3"),
#'     Relation_to_UCSC_CpG_Island = c("N_Shore","S_Shelf","N_Shore"),
#'     UCSC_CpG_Islands_Name = c("chrY:11-22","chrY:22-33","chrY:55-66"),
#'     UCSC_RefGene_Accession = c("NM_00116441;NR_00155","NM_00468","NR_00292"),
#'     Chromosome_36 = c("Y","Y","Y"),
#'     Coordinate_36 = c("9973356","21151183","20123684"),
#'     UCSC_RefGene_Name = c("ARHGEF10L;HIF3A","HIF3A","RNF17"),
#'     UCSC_RefGene_Group = c("Body;TSS1500","Body","Body"),
#'     stringsAsFactors=FALSE)
#'
#' Expressions <- data.frame(
#'     sample = c("ARHGEF10L", "HIF3A", "RNF17"),
#'     TCGA1 = c(-0.746592469762, -0.753826336325, -0.531035005853),
#'     TCGA2 = c(-0.363492469762, -2.67922633632, -0.531035005853),
#'     TCGA3 = c(-2.97279246976, 0.679073663675, -4.06262633632),
#'     TCGA4 = c(-0.613189285993,-0.545547600468, 0.59499523215),
#'     TCGA5 = c(-0.602787155677, -1.56558715568, 3.11398609442),
#'     TCGA6 = c(0.323339547405, 0.268878821693, -1.62560252829),
#'     stringsAsFactors=FALSE)
#'
#' Methylation <- data.frame(
#'     sample = c("cg1", "cg2", "cg3"),
#'     TCGA1 = c(0.1290, 0.2394, 0.0913),
#'     TCGA2 = c(0.9214, 0.0193, 0.0281),
#'     TCGA3 = c(0.5925, 0.8948, 0.0720),
#'     TCGA4 = c(0.9643, 0.7934, 0.0354),
#'     TCGA5 = c(0.1298, 0.0243, 0.8296),
#'     TCGA6 = c(0.8508, 0.8952, 0.9893),stringsAsFactors=FALSE)
#'
#' epimethex.analysis(Expressions, Annotations, Methylation, 1, 5, 2, FALSE)
#'
#' @source Functions
#'
#' @return files excel
#'
#' @export
epimethex.analysis <- function(dfPancan, dfGpl, dfMethylation,
        minRangeGene, maxRangeGene, numCores, dataLinear) {

    NUM_ROW_ADDED_FROM_ANALYSIS <- 9
    NUM_ROW_ADDED_FROM_ANALYSIS_ISLANDS_POSITIONSCG <- 17
    # dfPancan <- read.table(path_dfPancan, header = TRUE,
    #     sep = "\t")[minRangeGene:maxRangeGene,]
    dfPancan <- dfPancan[minRangeGene:maxRangeGene,]

    # folder name where xml files will be saved
    nameFD <- paste(as.character(minRangeGene), as.character(maxRangeGene),
        sep = "-")

    dir.create(file.path(getwd(), nameFD))

    # the first 37 rows of probset dataset are comments.
    # The following istruction reads only columns of interest
    # dfGpl <- data.table::fread(path_dfGpl,sep = "\t", header = TRUE,skip = 37,
    #     na.strings = c("", "NA"))[-1, c('ID', 'Relation_to_UCSC_CpG_Island',
    #     'UCSC_CpG_Islands_Name', 'UCSC_RefGene_Accession', 'Chromosome_36',
    #     'Coordinate_36', 'UCSC_RefGene_Name', 'UCSC_RefGene_Group' )]

    # create a new column `island` with the two columns collapsed together
    dfGpl$island <- apply(dfGpl[, c('Relation_to_UCSC_CpG_Island',
        'UCSC_CpG_Islands_Name')] ,1 , paste , collapse = "_")

    # remove the unnecessary columns
    dfGpl <-dfGpl[, -which(names(dfGpl) %in% c('Relation_to_UCSC_CpG_Island',
        'UCSC_CpG_Islands_Name'))]

    #dfMethylation <- read.table(path_dfMethylation, sep = "\t", header = TRUE)

    # transpose all but the first column (name)
    dfPancan2 <- as.data.frame(t(dfPancan[,-1]))
    colnames(dfPancan2) <- dfPancan$sample

    #remove all genes that have all values equal to zero
    dfPancan2 <- dfPancan2[, which(!apply(dfPancan2 == 0, 2, all))]

    #indexTCGA contains all the sorted patient's genes
    indexTCGA <- data.frame(sapply(seq_along(dfPancan2), function(x) {
        row.names(dfPancan2[order(dfPancan2[x], decreasing = TRUE),
        x, drop = FALSE])}))

    colnames(indexTCGA) <- colnames(dfPancan2)
    dfPancan2 <- apply(dfPancan2, 2, sort, decreasing = TRUE)

    #calculate quantili of all genes into dfPancan2
    if (requireNamespace("matrixStats", quietly = TRUE)) {
        quantili <- data.frame(matrixStats::colQuantiles(dfPancan2,
        probs = c(0.33, 0.66, 0.99)))
    }

    data.table::setnames(quantili, c("perc33", "perc66", "perc99"))
    quantili <- as.data.frame(t(quantili))

    dfPancan2 <- as.data.frame(dfPancan2)

    # [START] these istructions divide into equal parts the dfPancan2
    index <- which.max(apply(dfPancan2, 2, function(x) length(unique(x))))[[1]]

    dfPancan2$variable <- with( dfPancan2, RcmdrMisc::bin.var(
        dfPancan2[, index], bins = 3, method = 'proportions',
        labels = c('D', 'M', 'UP')))

    # [END]
    lengthDfpancan <- nrow(dfPancan2)

    # [START] calculate the means of all values of genes that are UP,
    # Medium, Down
    meanUP <-apply(dfPancan2[which(dfPancan2$variable %in% "UP"),
        -ncol(dfPancan2)], 2, mean)

    meanMID <- apply(dfPancan2[which(dfPancan2$variable %in% "M"),
        -ncol(dfPancan2)], 2, mean)

    meanDOWN <- apply(dfPancan2[which(dfPancan2$variable %in% "D"),
        -ncol(dfPancan2)], 2, mean)
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

    #Compute combinations fold change (up vs mid, up vs down, etc..)
    dfPancan2 <- calcFC(dfPancan2, dataLinear)

    dimDFpancan <- dim(dfPancan2[, -ncol(dfPancan2)])[2]

    # Create cluster with desired number of cores
    cl <- parallel::makeCluster(numCores)

    # Register cluster
    registerDoParallel(cl)

    indexTmp <- NULL
    #[START] calculate t-test for all genes
    resultUPvsMID <-foreach(indexTmp = seq_len(dimDFpancan), .combine = rbind)%dopar%{
        calculateTtest(dfPancan2[which(dfPancan2$variable %in% "UP"), indexTmp],
        dfPancan2[which(dfPancan2$variable %in% "M"), indexTmp], dataLinear)
    }

    resultUPvsDOWN <-foreach(indexTmp = seq_len(dimDFpancan), .combine = rbind)%dopar%{
        calculateTtest(dfPancan2[which(dfPancan2$variable %in% "UP"), indexTmp],
        dfPancan2[which(dfPancan2$variable %in% "D"), indexTmp], dataLinear)
    }

    resultMIDvsDOWN <-
        foreach(indexTmp = seq_len(dimDFpancan), .combine = rbind)%dopar%{
        calculateTtest(dfPancan2[which(dfPancan2$variable %in% "M"), indexTmp],
        dfPancan2[which(dfPancan2$variable %in% "D"), indexTmp], dataLinear)
    }

    parallel::stopCluster(cl)
    #[END]

    resultUPvsMID <- as.data.frame(t(resultUPvsMID))
    resultUPvsDOWN <- as.data.frame(t(resultUPvsDOWN))
    resultMIDvsDOWN <- as.data.frame(t(resultMIDvsDOWN))

    colnames(resultUPvsMID) <- colnames(dfPancan2[, -ncol(dfPancan2)])
    colnames(resultUPvsDOWN) <- colnames(dfPancan2[, -ncol(dfPancan2)])
    colnames(resultMIDvsDOWN) <- colnames(dfPancan2[, -ncol(dfPancan2)])

    dfPancan2 <- plyr::rbind.fill(dfPancan2,
        resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)

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

    rowNames <- c("meanDown", "meanMedium", "meanUP", "perc33", "perc66",
        "perc99", "fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN", "ttest_UPvsMID",
        "pvalue_UPvsMID", "ttest_UPvsDOWN", "pvalue_UPvsDOWN",
        "ttest_MIDvsDOWN", "pvalue_MIDvsDOWN" )

    j <- 1
    for (i in ROW_MEAN_DOWN:ROW_PVALUE_MIDvsDOWN) {
        row.names(dfPancan2)[i] <- rowNames[j]
        j <- j + 1
    }

    # [START] the dfAnnotations contain all cg with the appropriate genes,
    # islands, position of body and position into chromosome
    s <- strsplit(dfGpl$UCSC_RefGene_Name, split = ";")
    s1 <- strsplit(dfGpl$UCSC_RefGene_Accession, split = ";")
    s2 <- strsplit(dfGpl$UCSC_RefGene_Group, split = ";")

    dfAnnotations <- data.frame( cg = rep(dfGpl$ID, sapply(s, length)),
        Chromosome_36 = rep(dfGpl$Chromosome_36, sapply(s, length)),
        Coordinate_36 = rep(dfGpl$Coordinate_36, sapply(s, length)),
        gene = unlist(s), UCSC_RefGene_Accession = unlist(s1),
        position = unlist(s2), island = rep(dfGpl$island, sapply(s, length))
    )
    #[END]

    remove(s, s1, s2, dfGpl)
    gc()

    #remove genes and cg that aren't into dfPancan2 from dfAnnotations
    dfAnnotations <- subset(dfAnnotations,
        as.character(dfAnnotations$gene) %in% names(dfPancan2))

    #remove genes and cg that aren't into dfPancan2 from dfMethylation
    dfMethylation <-subset(dfMethylation,
        as.character(dfMethylation$sample) %in% dfAnnotations$cg)

    maxOccurence <- max(as.data.frame(table(unlist(dfAnnotations$gene)))$Freq)

    #reorder genes related CG
    dfAnnotations <-dfAnnotations %>% plyr::arrange(dfAnnotations$gene,
        dfAnnotations$Coordinate_36)

    #dfCGunique contains the unique couple (CG-position)
    dfCGunique <- dfAnnotations[!duplicated(dfAnnotations[, c(1, 6)]),]

    dfMethylation$sample <- as.factor(dfMethylation$sample)
    dfMethylation <- as.data.frame(dfMethylation)
    rownames(dfMethylation) <- dfMethylation$sample
    data.table::setkey(data.table::as.data.table(dfMethylation), sample)
    gc()

    tmp <- tapply(lapply(seq_len(nrow(dfCGunique)),function (i)
        (dfMethylation[as.vector(dfCGunique[i, 1]),
        as.vector(indexTCGA[,as.vector( dfCGunique[i, "gene"])])])),
        factor(dfCGunique[, "gene"]), function (x) { unname(unlist(x))})

    max.rows <- max(sapply(tmp, length))

    # DFCGorder is a dataframe that contains all CG values
    # ordered by indexTCGA dataframe
    DFCGorder <-do.call(cbind, lapply(tmp, function(x) {
        length(x) <- max.rows
        return (x)
    }))

    valExprGene <- dfPancan2[c(ROW_FC_UPvsMID:ROW_FC_MIDvsDOWN,
        ROW_PVALUE_UPvsMID, ROW_PVALUE_UPvsDOWN, ROW_PVALUE_MIDvsDOWN),
        -ncol(dfPancan2)]

    positions <- as.vector(unique(dfCGunique$position))
    genes <- colnames(DFCGorder)
    islands <- as.vector(unique(dfCGunique$island))
    #remove islands that have values "NA_NA"
    islands <- Filter(function(x) ! any(grepl("NA_NA", x)), islands)

    stratification <- as.data.frame(rep(c("UP", "Medium", "Down"),
        each = lengthDfpancan/3))
    colnames(stratification) <- "stratification"

    gc()
    lengthGens <- length(genes)
    lengthPositions <- length(positions)
    lengthIslands <- length(islands)

    # Create cluster with desired number of cores
    cl <- parallel::makeCluster(numCores)

    # Register cluster
    registerDoParallel(cl)

    parallel::clusterCall(cl, requireNamespace, package = "miscTools",
        quietly = TRUE)

    parallel::clusterCall(cl, requireNamespace, package = "psych",
        quietly = TRUE)

    parallel::clusterCall(cl, requireNamespace, package = "plyr",
        quietly = TRUE)

    # [START] here it is created a dataframe containing
    # means, beta differences, and pvalues of the stratifications (UP,MED,DOWN)
    # of all CG in all the postions of all genes.

    mFinaleCGglobali <- foreach(i = seq_len(lengthGens), .combine = cbind) %dopar% {

        flag <- FALSE

        # [START] create the tempMatrix that will contain all the CG values
        # of the gene "i"
        tempMatrix <- data.frame(matrix(DFCGorder[, i],nrow = lengthDfpancan,
            ncol = length(DFCGorder[, i]) / lengthDfpancan ))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete from tempMatrix all the columns that are not associated
        # to CG in gene "i"

        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")

        # columnNA contains all the NA set columns indexes
        columnNA <- which(sapply(tempMatrix, function(x) all(is.na(x))))

        tempMatrix <-
            tempMatrix[, colSums(is.na(tempMatrix)) != nrow(tempMatrix)]

        if (length(tempMatrix) != 0) {
            tempMatrix <- cbind(tempMatrix, stratification)
            tempMatrix <- Analysis(tempMatrix,dataLinear)
            tempMatrix$stratification <- NULL
        } else{
            flag <- TRUE
        }

        # insert in tempMatrix columns that were deleted in the same position

        if (length(columnNA) != 0) {
            for (j in seq_len(length(columnNA))) {
                tempMatrix <-
                insertCol(as.matrix(tempMatrix), columnNA[[j]], v = NA)
            }
        }

        if (flag) {
            df <- as.data.frame(matrix(NA, nrow = 9, ncol = dim(tempMatrix)[2]))
            tempMatrix <- rbind(tempMatrix, df)
        }

        # [START] compute the correlation among gene expression data
        # and methylation data of associated CG

        DataexpressionGeneTmp <-
            as.data.frame(dfPancan2[c(seq_len(lengthDfpancan)), genes[i]])

        DataCG_Tmp <- as.data.frame(tempMatrix[c(seq_len(lengthDfpancan)), ])

        # compute the correlation test among gene expression data and CG data

        resultCorrTest <-
            psych::corr.test(DataexpressionGeneTmp, DataCG_Tmp, adjust = "none")

        tempMatrix <- rbind(tempMatrix,as.numeric(resultCorrTest$r),
            as.numeric(resultCorrTest$p))

        # [START] delete from tempMatrix the CG expression values, keeping
        # only the followin rows: "medianDown", "medianMedium", "medianUP",
        # "bd_UPvsMID", "bd_UPvsDOWN", "bd_MIDvsDOWN", "pvalue_UPvsMID",
        # "pvalue_UPvsDOWN", "pvalue_MIDvsDOWN", "pearson_correlation",
        #  "pvalue_pearson_correlation"

        if (dim(tempMatrix)[1] > lengthDfpancan) {
            tempMatrix <- as.data.frame(tempMatrix[-c(seq_len(lengthDfpancan)), ])
        }
        #[END]

        # insert into temp matrix "tempMatrix" the following values of the
        # gene[i]: "fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN",
        # "pvalue_UPvsMID", "pvalue_UPvsDOWN", "pvalue_MIDvsDOWN"

        tempMatrix <- data.frame(sapply(tempMatrix, c,
            unlist(valExprGene[, genes[i]])), row.names = NULL)

        # [START] change the columns name associating the CG name with the
        # appropriate position gene (for example cg27394127_TSS1500_BFAR)

        colnames(tempMatrix) <-
            paste(as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]),
            "cg"]),
            as.character(dfCGunique[
                which(dfCGunique$gene %in% genes[i]), "position"]),
            genes[i], sep = "_")
        #[END]

        remove(DataexpressionGeneTmp, DataCG_Tmp, resultCorrTest)

        #print the matrix1 into mFinaleCGglobali
        tempMatrix
    }

    #[END]

    mFinaleCGglobali <- setRowNames(mFinaleCGglobali)

    # [START] extract the CG names, then recover the CG associated isles
    tmp <- strsplit(colnames(mFinaleCGglobali), split = "_")
    tmp <- unlist(lapply(tmp, '[[', 1))

    mFinaleCGglobali <- rbind(mFinaleCGglobali,
        as.character(dfCGunique[
        which(dfCGunique$cg %in% as.array(tmp)), "island"]))

    row.names(mFinaleCGglobali)[nrow(mFinaleCGglobali)] <- "island"
    #[END]
    mFinaleCGglobali <- t(mFinaleCGglobali)

    xlsx::write.xlsx( mFinaleCGglobali,
        file.path(getwd(), nameFD, "CG_data_individually.xlsx"),
        sheetName = "Sheet1")

    remove(mFinaleCGglobali, tmp)
    gc()

    mFinaleCGunificati <- foreach(i = seq_len(lengthGens), .combine = cbind) %dopar% {

        # [START] create tempMatrix containing the CG
        #values associated to gene[i]
        tempMatrix <- data.frame(matrix( DFCGorder[, i], nrow = lengthDfpancan,
            ncol = length(DFCGorder[, i]) / lengthDfpancan))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete in tempMatrix the columns
        # not associated in any CG of gene[i]
        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")

        tempMatrix <-
            as.data.frame(tempMatrix[,
            colSums(is.na(tempMatrix))!=nrow(tempMatrix)])

        num_CG <- length(tempMatrix)

        if (dim(tempMatrix)[2] != 0) {
            if (dim(tempMatrix)[2] != 1) {
                tempMatrix <- stack(tempMatrix)
                tempMatrix <- as.matrix(tempMatrix[,-2])
            }
            num_row_m4 <- nrow(tempMatrix)
            tempMatrix <- cbind(tempMatrix, stratification)
            colnames(tempMatrix) <- c("value", "stratification")
            tempMatrix <- Analysis(tempMatrix,dataLinear)
            tempMatrix <- as.data.frame(tempMatrix[,-2])

            names(tempMatrix) <- "value"

            # [START] compute the correlation among gene[i] expression data and
            # methylation data of associated CG
            dfTmp <-
                as.data.frame(rep(dfPancan2[c(seq_len(lengthDfpancan)), genes[i]],
                num_CG))

            m4Tmp <- as.data.frame(tempMatrix[c(seq_len(num_row_m4)), ])
            resultCorrTest <- psych::corr.test(dfTmp, m4Tmp, adjust = "none")

            tempMatrix <- rbind(tempMatrix, as.numeric(resultCorrTest$r),
                as.numeric(resultCorrTest$p))
            #[END]

            tempMatrix <- as.data.frame(tempMatrix[-c(seq_len(num_row_m4)),])
            tempMatrix <- data.frame(sapply(tempMatrix, c,
                unlist(valExprGene[, genes[i]])), row.names = NULL)
            colnames(tempMatrix) <- paste("CG", genes[i], sep = "_")
            remove(dfTmp, m4Tmp, resultCorrTest)

            tempMatrix
        } else{
            tempMatrix <- as.data.frame(matrix(NA,
            nrow = NUM_ROW_ADDED_FROM_ANALYSIS_ISLANDS_POSITIONSCG, ncol = 1))
            colnames(tempMatrix) <- paste("CG", genes[i], sep = "_")

            tempMatrix
        }
    }

    mFinaleCGunificati <- setRowNames(mFinaleCGunificati)
    mFinaleCGunificati <- t(mFinaleCGunificati)

    xlsx::write.xlsx( mFinaleCGunificati,
        file.path(getwd(), nameFD, "CG_of_a_genes.xlsx"),
            sheetName = "Sheet1")

    ###CG group by gene position

    mFinaleCGposition <- foreach(i = seq_len(lengthGens), .combine = cbind) %dopar% {

        # [START] create tempMatrix containing the CG values
        # associated to gene[i]
        tempMatrix <- data.frame(matrix( DFCGorder[, i], nrow = lengthDfpancan,
            ncol = length(DFCGorder[, i]) / lengthDfpancan ))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete from tempMatrix the columns not associated
        # to any CG in gene[i]

        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")

        AnalysisIslands_PositionsCG(lengthPositions,i,positions,"position",
            dfCGunique,genes,tempMatrix,stratification,dfPancan2,valExprGene,
            dataLinear,lengthDfpancan )
    }

    if(plyr::empty(mFinaleCGposition) == FALSE){

        mFinaleCGposition <- setRowNames(mFinaleCGposition)
        mFinaleCGposition <- t(mFinaleCGposition)
        xlsx::write.xlsx( mFinaleCGposition, file.path(getwd(), nameFD,
            "CG_by_position.xlsx"), sheetName = "Sheet1")
    }
    else{
        warning("CG_by_position.xlsx is not generated")
    }

    #CG group by island CpG
    mFinaleCGisland <- foreach(i = seq_len(lengthGens), .combine = cbind) %dopar% {

        # [START] create tempMatrix that contains the
        # CG values associated to gene[i]
        tempMatrix <- data.frame(matrix(DFCGorder[, i], nrow = lengthDfpancan,
            ncol = length(DFCGorder[, i]) / lengthDfpancan ))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete from tempMatrix the columns
        # not associated to any CG of gene[i]
        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "_")

        AnalysisIslands_PositionsCG(lengthIslands,i,islands,"island",dfCGunique,
            genes,tempMatrix,stratification,dfPancan2,valExprGene,
            dataLinear,lengthDfpancan)
    }

    if(plyr::empty(mFinaleCGisland) == FALSE) {

        mFinaleCGisland <- setRowNames(mFinaleCGisland)
        mFinaleCGisland <- t(mFinaleCGisland)
        xlsx::write.xlsx(mFinaleCGisland, file.path(getwd(), nameFD,
            "CG_by_islands.xlsx"), sheetName = "Sheet1")
    }
    else{
        warning("CG_by_islands.xlsx is not generated")
    }
}
