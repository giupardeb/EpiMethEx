#' Create a Epimethex function
#'
#' @param dfExpressions data frame, expression of data
#' @param dfGpl data frame, Annotations of data
#' @param dfMethylation data frame, Methylation of data
#' @param minRangeGene Numeric, lowerbound of genes
#' @param maxRangeGene Numeric, upperbound of genes
#' @param numCores Numeric, is the number of cores that you can use
#' @param dataGenesLinear Logic, determines if genetic data are linear
#' @param testExpression Logic, determines the test to apply on expression
#' dataset. If TRUE will apply t-student test,
#' otherwise will apply Kolmogorov-Smirnov test
#' @param testMethylation Logic, determines the test to apply on methylation
#' dataset. If TRUE will apply t-student test,
#' otherwise will apply Kolmogorov-Smirnov test
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
#' @importFrom MultiAssayExperiment assay
#'
#' @examples
#'
#' Annotations <- data.frame(
#' ID = c("cg11663302","cg01552731", "cg09081385"),
#' Relation_to_UCSC_CpG_Island = c("Island","N_Shore","N_Shore"),
#' UCSC_CpG_Islands_Name = c("chr1:18023481-18023792","chr19:46806998-46807617",
#' "chr12:120972167-120972447"),
#' UCSC_RefGene_Accession = c("NM_001011722","NM_152794","NM_014868"),
#' Chromosome_36 = c("1","19","12"),
#' Coordinate_36 = c("17896255","51498747","119456453"),
#' UCSC_RefGene_Name = c("ARHGEF10L","HIF3A","RNF10"),
#' UCSC_RefGene_Group =c("Body","1stExon","TSS200"),
#' stringsAsFactors=FALSE)
#'
#' Expressions <- data.frame(
#' 'sample' = c("ARHGEF10L", "HIF3A", "RNF10"),
#' 'TCGA-YD-A89C-06' = c(-0.746592469762, -0.753826336325, 0.4953280),
#' 'TCGA-Z2-AA3V-06' = c(0.578807530238, -2.30662633632, 0.1023280),
#' 'TCGA-EB-A3Y6-01' = c(-0.363492469762, -2.67922633632, -0.6147720),
#' 'TCGA-EE-A3JA-06' = c(-2.97279246976, -3.61932633632, 0.02932801),
#' 'TCGA-D9-A4Z2-01' = c(-0.128492469762, 0.679073663675, 0.4017280),
#' 'TCGA-D3-A51G-06' = c(-0.4299925, -4.0626263, -1.0136720),
#' stringsAsFactors=FALSE)
#'
#' Methylation <- data.frame(
#' 'sample' = c("cg11663302", "cg01552731", "cg09081385"),
#' 'TCGA-YD-A89C-06' = c(0.9856, 0.7681, 0.0407),
#' 'TCGA-Z2-AA3V-06' = c(0.9863, 0.8551, 0.0244),
#' 'TCGA-EB-A3Y6-01' = c(0.9876, 0.6473, 0.028),
#' 'TCGA-EE-A3JA-06' = c(0.9826, 0.4587, 0.0343),
#' 'TCGA-D9-A4Z2-01' = c(0.9881, 0.8509, 0.0215),
#' 'TCGA-D3-A51G-06' = c(0.9774, 0.813, 0.0332),
#' stringsAsFactors=FALSE)
#'
#' epimethex.analysis(Expressions, Annotations, Methylation, 1, 3, 2,
#' TRUE, TRUE, FALSE)
#'
#' @source Functions
#'
#' @return csv files
#'
#' @export
epimethex.analysis <- function(dfExpressions, dfGpl, dfMethylation,
        minRangeGene, maxRangeGene, numCores, dataGenesLinear,
        testExpression, testMethylation) {

    NUM_ROW_ADDED_FROM_ANALYSIS <- 9
    NUM_ROW_ADDED_FROM_ANALYSIS_ISLANDS_POSITIONSCG <- 17


    if(typeof(dfExpressions) == "S4"){
        tmp1<-assay(dfExpressions)
        rm(dfExpressions)
        gc()
        mode(tmp1) = "numeric"
        dfExpressions <- data.frame(tmp1)
    }

    dfExpressions <- dfExpressions[minRangeGene:maxRangeGene,]
    # folder name where xml files will be saved
    nameFD <- paste(as.character(minRangeGene), as.character(maxRangeGene),
        sep = "-")

    dir.create(file.path(getwd(), nameFD),showWarnings = FALSE)

    dfGpl<-dfGpl[!(is.na(dfGpl$UCSC_RefGene_Name)|dfGpl$UCSC_RefGene_Name==""),]
    # create a new column `island` with the two columns collapsed together
    dfGpl$island <- apply(dfGpl[, c('Relation_to_UCSC_CpG_Island',
        'UCSC_CpG_Islands_Name')] ,1 , paste , collapse = "_")

    # remove the unnecessary columns
    dfGpl <-dfGpl[, -which(names(dfGpl) %in% c('Relation_to_UCSC_CpG_Island',
        'UCSC_CpG_Islands_Name'))]

    # transpose all but the first column (name)
    dfExpressions2 <- as.data.frame(t(dfExpressions[,-1]))
    if("sample" %in% colnames(dfExpressions)) {
        colnames(dfExpressions2) <- dfExpressions$sample
    }else{
        colnames(dfExpressions2) <- rownames(dfExpressions)
    }


    #remove all genes that have all values equal to zero
    dfExpressions2 <- dfExpressions2[, which(!apply(
        dfExpressions2 == 0, 2, all))]
    dfExpressions2<-na.omit(dfExpressions2)

    #indexTCGA contains all the sorted patient's genes
    indexTCGA <- data.frame(sapply(seq_along(dfExpressions2), function(x) {
        row.names(dfExpressions2[order(dfExpressions2[x], decreasing = TRUE),
        x, drop = FALSE])}))

    colnames(indexTCGA) <- colnames(dfExpressions2)
    dfExpressions2 <- apply(dfExpressions2, 2, sort, decreasing = TRUE)

    #calculate quantili of all genes into dfExpressions2
    quantili <- data.frame(matrixStats::colQuantiles(dfExpressions2,
        probs = c(0.33, 0.66, 0.99)))

    data.table::setnames(quantili, c("perc33", "perc66", "perc99"))
    quantili <- as.data.frame(t(quantili))

    dfExpressions2 <- as.data.frame(dfExpressions2)

    # [START] these istructions divide into equal parts the dfExpressions2
    index <- which.max(apply(dfExpressions2, 2,
        function(x) length(unique(x))))[[1]]

    dfExpressions2$variable <- with( dfExpressions2, RcmdrMisc::bin.var(
        dfExpressions2[, index], bins = 3, method = 'proportions',
        labels = c('D', 'M', 'UP')))

    # [END]
    lengthdfExpressions <- nrow(dfExpressions2)

    # [START] calculate the means of all values of genes that are UP,
    # Medium, Down
    meanUP <-apply(dfExpressions2[which(dfExpressions2$variable %in% "UP"),
        -ncol(dfExpressions2)], 2, mean)

    meanMID <- apply(dfExpressions2[which(dfExpressions2$variable %in% "M"),
        -ncol(dfExpressions2)], 2, mean)

    meanDOWN <- apply(dfExpressions2[which(dfExpressions2$variable %in% "D"),
        -ncol(dfExpressions2)], 2, mean)
    # [END]
    dfExpressions2 <- rbind(dfExpressions2, meanDOWN, meanMID, meanUP)

    remove(meanUP, meanMID, meanDOWN)

    dfExpressions2 <- plyr::rbind.fill(dfExpressions2, quantili)

    remove(quantili)

    ROW_MEAN_DOWN <- nrow(dfExpressions2) - 5
    ROW_MEAN_MEDIUM <- nrow(dfExpressions2) - 4
    ROW_MEAN_UP <- nrow(dfExpressions2) - 3
    ROW_PERC_33 <- nrow(dfExpressions2) - 2
    ROW_PERC_66 <- nrow(dfExpressions2) - 1
    ROW_PERC_99 <- nrow(dfExpressions2)
    gc()

    #Compute combinations fold change (up vs mid, up vs down, etc..)
    dfExpressions2 <- calcFC(dfExpressions2, dataGenesLinear)

    dimdfExpressions <- dim(dfExpressions2[, -ncol(dfExpressions2)])[2]

    # Create cluster with desired number of cores
    cl <- parallel::makeCluster(numCores)

    # Register cluster
    registerDoParallel(cl)

    indexTmp <- NULL
    #[START] calculate t-test for all genes
    resultUPvsMID <-
    foreach(indexTmp = seq_len(dimdfExpressions), .combine = rbind)%dopar%{
        calculateTtest(dfExpressions2[which(dfExpressions2$variable %in% "UP"),
            indexTmp],
        dfExpressions2[which(dfExpressions2$variable %in% "M"), indexTmp],
        testExpression)
    }

    resultUPvsDOWN <-
    foreach(indexTmp = seq_len(dimdfExpressions), .combine = rbind)%dopar%{
        calculateTtest(dfExpressions2[which(dfExpressions2$variable %in% "UP"),
            indexTmp],
        dfExpressions2[which(dfExpressions2$variable %in% "D"), indexTmp],
        testExpression)
    }

    resultMIDvsDOWN <-
    foreach(indexTmp = seq_len(dimdfExpressions), .combine = rbind)%dopar%{
        calculateTtest(dfExpressions2[which(dfExpressions2$variable %in% "M"),
            indexTmp],
        dfExpressions2[which(dfExpressions2$variable %in% "D"), indexTmp],
        testExpression)
    }

    parallel::stopCluster(cl)
    #[END]

    resultUPvsMID <- as.data.frame(t(resultUPvsMID))
    resultUPvsDOWN <- as.data.frame(t(resultUPvsDOWN))
    resultMIDvsDOWN <- as.data.frame(t(resultMIDvsDOWN))

    colnames(resultUPvsMID) <-colnames(dfExpressions2[, -ncol(dfExpressions2)])
    colnames(resultUPvsDOWN) <-colnames(dfExpressions2[,-ncol(dfExpressions2)])
    colnames(resultMIDvsDOWN) <-colnames(dfExpressions2[,-ncol(dfExpressions2)])

    dfExpressions2 <- plyr::rbind.fill(dfExpressions2,
        resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)

    remove(resultUPvsMID, resultUPvsDOWN, resultMIDvsDOWN)
    gc()

    ROW_FC_UPvsMID <- nrow(dfExpressions2) - 8
    ROW_FC_UPvsDOWN <- nrow(dfExpressions2) - 7
    ROW_FC_MIDvsDOWN <- nrow(dfExpressions2) - 6
    ROW_TTEST_UPvsMID <- nrow(dfExpressions2) - 5
    ROW_PVALUE_UPvsMID <- nrow(dfExpressions2) - 4
    ROW_TTEST_UPvsDOWN <- nrow(dfExpressions2) - 3
    ROW_PVALUE_UPvsDOWN <- nrow(dfExpressions2) - 2
    ROW_TTEST_MIDvsDOWN <- nrow(dfExpressions2) - 1
    ROW_PVALUE_MIDvsDOWN <- nrow(dfExpressions2)

    rowNames <- c("meanDown", "meanMedium", "meanUP", "perc33", "perc66",
        "perc99", "fc_UPvsMID", "fc_UPvsDOWN", "fc_MIDvsDOWN", "ttest_UPvsMID",
        "pvalue_UPvsMID", "ttest_UPvsDOWN", "pvalue_UPvsDOWN",
        "ttest_MIDvsDOWN", "pvalue_MIDvsDOWN" )

    j <- 1
    for (i in ROW_MEAN_DOWN:ROW_PVALUE_MIDvsDOWN) {
        row.names(dfExpressions2)[i] <- rowNames[j]
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

    #remove genes and cg that aren't into dfExpressions2 from dfAnnotations
    dfAnnotations <- subset(dfAnnotations,
        as.character(dfAnnotations$gene) %in% names(dfExpressions2))

    if(typeof(dfMethylation) == "S4"){

        dftmp1<-MultiAssayExperiment::assay(dfMethylation)
        mode(dftmp1) = "numeric"
        rm(dfMethylation)
        dfMethylation <- as.data.frame(dftmp1)
        rm(dftmp1)
        dfMethylation <- cbind(sample=rownames(dfMethylation),
            dfMethylation, row.names = NULL)
    }

    dfMethylation<-na.omit(dfMethylation)
    #remove genes and cg that aren't into dfExpressions2 from dfMethylation
    colnames(dfMethylation)[1] <- "sample"
    dfMethylation <-subset(dfMethylation,
        as.character(dfMethylation$sample) %in% as.character(dfAnnotations$cg))

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

    valExprGene <- dfExpressions2[c(ROW_FC_UPvsMID:ROW_FC_MIDvsDOWN,
        ROW_PVALUE_UPvsMID, ROW_PVALUE_UPvsDOWN, ROW_PVALUE_MIDvsDOWN),
        -ncol(dfExpressions2)]

    positions <- as.vector(unique(dfCGunique$position))
    genes <- colnames(DFCGorder)
    islands <- as.vector(unique(dfCGunique$island))
    #remove islands that have values "NA_NA" or "_" value
    islands <- Filter(function(x) ! any(grepl("NA_NA|^_", x)), islands)

    stratification <- as.data.frame(rep(c("UP", "Medium", "Down"),
        each = lengthdfExpressions/3))
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

    mFinaleCGglobali <-
    foreach(i = seq_len(lengthGens), .combine = cbind)%dopar% {

        flag <- FALSE

        # [START] create the tempMatrix that will contain all the CG values
        # of the gene "i"
        tempMatrix <- data.frame(matrix(DFCGorder[, i],nrow=lengthdfExpressions,
            ncol = length(DFCGorder[, i]) / lengthdfExpressions ))

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
            tempMatrix <- Analysis(tempMatrix,testMethylation)
            tempMatrix$stratification <- NULL
        } else{
            flag <- TRUE
        }

        # insert in tempMatrix columns that were deleted in the same position

        if (length(columnNA) != 0) {
            for (j in seq_len(length(columnNA))) {
                tempMatrix<-insertCol(as.matrix(tempMatrix),columnNA[[j]],v=NA)
            }
        }

        if (flag) {
            df <- as.data.frame(matrix(NA, nrow = 9, ncol = dim(tempMatrix)[2]))
            tempMatrix <- rbind(tempMatrix, df)
        }

        # [START] compute the correlation among gene expression data
        # and methylation data of associated CG

        DataexpressionGeneTmp <-
            as.data.frame(dfExpressions2[c(seq_len(lengthdfExpressions)),
                genes[i]])

        DataCG_Tmp <-
            as.data.frame(tempMatrix[c(seq_len(lengthdfExpressions)), ])

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

        if (dim(tempMatrix)[1] > lengthdfExpressions) {
            tempMatrix <-
                as.data.frame(tempMatrix[-c(seq_len(lengthdfExpressions)),])
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
    write.csv(mFinaleCGglobali,file.path(getwd(),
        nameFD,"CG_data_individually.csv"))

    remove(mFinaleCGglobali, tmp)
    gc()

    mFinaleCGunificati <-
    foreach(i = seq_len(lengthGens), .combine = cbind) %dopar% {

        # [START] create tempMatrix containing the CG
        #values associated to gene[i]
        tempMatrix<-data.frame(matrix( DFCGorder[, i], nrow=lengthdfExpressions,
            ncol = length(DFCGorder[, i]) / lengthdfExpressions))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete in tempMatrix the columns
        # not associated in any CG of gene[i]
        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "~")

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
            tempMatrix <- Analysis(tempMatrix,testMethylation)
            tempMatrix <- as.data.frame(tempMatrix[,-2])

            names(tempMatrix) <- "value"

            # [START] compute the correlation among gene[i] expression data and
            # methylation data of associated CG
            dfTmp <- as.data.frame(rep(dfExpressions2[
                c(seq_len(lengthdfExpressions)), genes[i]], num_CG))

            m4Tmp <- as.data.frame(tempMatrix[c(seq_len(num_row_m4)), ])
            resultCorrTest <- psych::corr.test(dfTmp, m4Tmp, adjust = "none")

            tempMatrix <- rbind(tempMatrix, as.numeric(resultCorrTest$r),
                as.numeric(resultCorrTest$p))
            #[END]

            tempMatrix <- as.data.frame(tempMatrix[-c(seq_len(num_row_m4)),])
            tempMatrix <- data.frame(sapply(tempMatrix, c,
                unlist(valExprGene[, genes[i]])), row.names = NULL)
            colnames(tempMatrix) <- paste("CG", genes[i], sep = "~")
            remove(dfTmp, m4Tmp, resultCorrTest)

            tempMatrix
        } else{
            tempMatrix <- as.data.frame(matrix(NA,
            nrow = NUM_ROW_ADDED_FROM_ANALYSIS_ISLANDS_POSITIONSCG, ncol = 1))
            colnames(tempMatrix) <- paste("CG", genes[i], sep = "~")

            tempMatrix
        }
    }

    mFinaleCGunificati <- setRowNames(mFinaleCGunificati)
    mFinaleCGunificati <- t(mFinaleCGunificati)
    write.csv(mFinaleCGunificati,file.path(getwd(), nameFD,"CG_of_a_genes.csv"))

    ###CG group by gene position

    mFinaleCGposition <-
    foreach(i = seq_len(lengthGens), .combine = cbind)%dopar% {

        # [START] create tempMatrix containing the CG values
        # associated to gene[i]
        tempMatrix<-data.frame(matrix( DFCGorder[, i], nrow=lengthdfExpressions,
            ncol = length(DFCGorder[, i]) / lengthdfExpressions ))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete from tempMatrix the columns not associated
        # to any CG in gene[i]

        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "~")

        AnalysisIslands_PositionsCG(lengthPositions,i,positions,"position",
            dfCGunique,genes,tempMatrix,stratification,dfExpressions2,
            valExprGene, testMethylation,lengthdfExpressions )
    }

    if(plyr::empty(mFinaleCGposition) == FALSE){

        mFinaleCGposition <- setRowNames(mFinaleCGposition)
        mFinaleCGposition <- t(mFinaleCGposition)
        write.csv(mFinaleCGposition,file.path(getwd(),
            nameFD,"CG_by_position.csv"))
    }else{
        warning("CG_by_position.csv is not generated")
    }

    #CG group by island CpG
    mFinaleCGisland <-
    foreach(i = seq_len(lengthGens), .combine = cbind)%dopar% {

        # [START] create tempMatrix that contains the
        # CG values associated to gene[i]
        tempMatrix<-data.frame(matrix(DFCGorder[, i], nrow=lengthdfExpressions,
            ncol = length(DFCGorder[, i]) / lengthdfExpressions ))

        colnames(tempMatrix) <-
            as.character(dfCGunique[which(dfCGunique$gene %in% genes[i]), "cg"])
        #[END]

        # [START] delete from tempMatrix the columns
        # not associated to any CG of gene[i]
        keep.cols <- names(tempMatrix) %in% NA
        tempMatrix <- tempMatrix [!keep.cols]
        #[END]

        colnames(tempMatrix) <- paste(names(tempMatrix), genes[i], sep = "~")

        AnalysisIslands_PositionsCG(lengthIslands,i,islands,"island",dfCGunique,
            genes,tempMatrix,stratification,dfExpressions2,valExprGene,
            testMethylation,lengthdfExpressions)
    }

    if(plyr::empty(mFinaleCGisland) == FALSE) {

        mFinaleCGisland <- setRowNames(mFinaleCGisland)
        mFinaleCGisland <- t(mFinaleCGisland)
        write.csv(mFinaleCGisland,file.path(getwd(),
            nameFD,"CG_by_islands.csv"))
    }
    else{
        warning("CG_by_islands.csv is not generated")
    }
}
