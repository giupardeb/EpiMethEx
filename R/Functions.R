#[START] this function calculates the Beta Difference
calcBetaDifference <- function(matrix3) {
    medianUP <- matrix3[nrow(matrix3),-ncol(matrix3)]
    medianMID <- matrix3[nrow(matrix3) - 1,-ncol(matrix3)]
    medianDOWN <- matrix3[nrow(matrix3) - 2,-ncol(matrix3)]

    if (dim(matrix3)[2] == 2) {
        bd_UPvsMID <- data.frame(value = medianUP - medianMID)
        bd_UPvsDOWN <- data.frame(value = medianUP - medianDOWN)
        bd_MIDvsDOWN <- data.frame(value = medianMID - medianDOWN)
    } else{
        bd_UPvsMID <- data.frame(medianUP - medianMID)
        bd_UPvsDOWN <- data.frame(medianUP - medianDOWN)
        bd_MIDvsDOWN <- data.frame(medianMID - medianDOWN)
        colnames(bd_UPvsMID) <- colnames(matrix3[,-ncol(matrix3)])
        colnames(bd_UPvsDOWN) <- colnames(matrix3[,-ncol(matrix3)])
        colnames(bd_MIDvsDOWN) <- colnames(matrix3[,-ncol(matrix3)])
    }

    matrix3 <- plyr::rbind.fill(matrix3, bd_UPvsMID, bd_UPvsDOWN, bd_MIDvsDOWN)
    return(matrix3)
}
#[END]

#[START] this function calculates the Fold Change
calcFC <- function(matrix2, flag) {
    df <- data.frame()

    meanUp <- matrix2[nrow(matrix2),-ncol(matrix2)]
    meanMid <- matrix2[nrow(matrix2) - 1,-ncol(matrix2)]
    meanDown <- matrix2[nrow(matrix2) - 2,-ncol(matrix2)]

    #flag = True, data are linear, else data are logarithmic
    if (flag) {
        fc_UPvsMID <- calculateLinearFC(meanUp, meanMid)
        fc_UPvsDOWN <-  calculateLinearFC(meanUp, meanDown)
        fc_MIDvsDOWN <- calculateLinearFC(meanMid, meanDown)
    } else{
        fc_UPvsMID <- calculateLogFC(meanUp, meanMid)
        fc_UPvsDOWN <- calculateLogFC(meanUp, meanDown)
        fc_MIDvsDOWN <-  calculateLogFC(meanMid, meanDown)
    }

    df <- rbind(df, fc_UPvsMID, fc_UPvsDOWN, fc_MIDvsDOWN)
    colnames(df) <- colnames(matrix2[, -ncol(matrix2)])

    matrix2 <- plyr::rbind.fill(matrix2, df)
    remove(df)

    return(matrix2)
}
#[END]

checkSign <- function(a, b) {
    return (sign(a) == sign(b))
}

#[START] this function calculate the Fold Change when data are logarithmic
calculateLogFC <- function(meanFirstGroup, meanSecondGroup) {
    fold <- 2 ^ (abs(meanFirstGroup - meanSecondGroup))
    fc <- ifelse(meanSecondGroup < meanFirstGroup, fold,-fold)
    return(fc)
}
#[END]

#[START] this function calculates the Fold Change when data are representd in linear way
calculateLinearFC <- function(meanFirstGroup, meanSecondGroup) {
    maxabs <- mapply(max, abs(meanFirstGroup), abs(meanSecondGroup))
    minabs <- mapply(min, abs(meanFirstGroup), abs(meanSecondGroup))
    max <- mapply(max, meanFirstGroup, meanSecondGroup)
    min <- mapply(min, meanFirstGroup, meanSecondGroup)
    Y <- checkSign(meanFirstGroup, meanSecondGroup)

    fc <- ifelse(Y == TRUE, maxabs / minabs, max - min)
    fold <- ifelse(meanSecondGroup < meanFirstGroup, fc,-fc)

    return(fold)
}
#[END]

# [START]
# This function analyzes global CG, computes the median, beta-difference and
# Kolmogorov-Smirnov test for the stratifications
Analysis <- function(matrix1) {
    medianUP <- apply(as.data.frame(
        matrix1[which(matrix1$stratification %in% "UP"), -ncol(matrix1)]), 2,
        median, na.rm = TRUE)

    medianMID <-apply(as.data.frame(
        matrix1[which(matrix1$stratification %in% "Medium"),-ncol(matrix1)]), 2,
        median, na.rm = TRUE)

    medianDOWN <-  apply(as.data.frame(
        matrix1[which(matrix1$stratification %in% "Down"),-ncol(matrix1)]), 2,
        median, na.rm = TRUE)

    if (dim(matrix1)[2] == 2) {
        colnames(matrix1) <- c("value", "stratification")
        colnamesDfTtest <- "value"
        dimM <- 1
    } else{
        dimM <- dim(matrix1[,-ncol(matrix1)])[2]
        colnamesDfTtest <- colnames(matrix1[, -ncol(matrix1)])
    }

    dfTtest <- data.table::setnames(data.frame(matrix(nrow = 3, ncol = dimM)),
        colnamesDfTtest)

    matrix1 <- rbind(matrix1, medianDOWN, medianMID, medianUP)
    remove(medianUP, medianMID, medianDOWN)

    matrix1 <- calcBetaDifference(matrix1)

    for (k in 1:dimM) {
        dfTtest[1, k] <- calculateTtest(matrix1[
            which(matrix1$stratification %in% "UP"), k],
            matrix1[which(matrix1$stratification %in% "Medium"), k], FALSE)

        dfTtest[2, k] <-calculateTtest(matrix1[
            which(matrix1$stratification %in% "UP"), k],
            matrix1[which(matrix1$stratification %in% "Down"), k], FALSE)

        dfTtest[3, k] <- calculateTtest(matrix1[
            which(matrix1$stratification %in% "Medium"), k],
            matrix1[which(matrix1$stratification %in% "Down"), k], FALSE)
    }

    matrix1 <- plyr::rbind.fill(matrix1, dfTtest)
    remove(dfTtest)
    return(matrix1)
}
# [END]

#Used for islands and positions cg grouping
AnalysisIslands_PositionsCG <- function(leng,index,position,column,dfCGunique,
    genes,tempMatrix,stratification,dfPancan2,valExprGene) {

    mFinale <- data.frame(matrix())

    for (k in 1:leng) {
        cg <- as.vector(dfCGunique[which(dfCGunique$gene %in% genes[index] &
            dfCGunique[, column] %in% position[k]), 1])

        if (length(cg) != 0) {
            cg <- paste(cg, genes[index], sep = "_")
            tempMatrix2 <- as.data.frame(tempMatrix[, cg])
            #remove all coloumns that have all values NA
            tempMatrix2 <- as.data.frame(tempMatrix2[,colSums(
                is.na(tempMatrix2)) != nrow(tempMatrix2)])

            if (dim(tempMatrix2)[2] > 1) {
                tempMatrix2 <- stack(tempMatrix2)
                tempMatrix2 <- as.matrix(tempMatrix2[, -2])
                num_row <- nrow(tempMatrix2)
                tempMatrix2 <- cbind(tempMatrix2, stratification)
                colnames(tempMatrix2) <- c("value", "stratification")

                tempMatrix2 <- Analysis(tempMatrix2)
                tempMatrix2 <- as.data.frame(tempMatrix2[, -2])

                # [START] computes the correlation among gene[i] expression data
                # and methylation data with associated CG

                mTmp <- as.data.frame(rep(dfPancan2[c(1:473), genes[index]],
                    length(cg)))

                mTmp1 <- as.data.frame(tempMatrix2[1:dim(mTmp)[1],])
                resultCorrTest <- psych::corr.test(mTmp, mTmp1, adjust = "none")

                #add pearson correlation and p-value
                tempMatrix2 <- rbind(tempMatrix2,as.numeric(resultCorrTest$r),
                    as.numeric(resultCorrTest$p))
                #[END]

                tempMatrix2 <- data.frame(sapply(tempMatrix2, c,
                    unlist(valExprGene[, genes[index]])), row.names = NULL)

                tempMatrix2 <- as.data.frame(tempMatrix2[-c(1:num_row), ])
                colnames(tempMatrix2) <- paste(position[k], genes[index],
                    sep = "_")
                mFinale <- cbind(mFinale, tempMatrix2)
            }
        }
    }
    return(subset(mFinale, select = -c(1)))
}

# [START] This function computes the t-student test in genes
# analysis when flag = F, and calculate Kolmogorov-Smirnov Tests otherwhise
calculateTtest <- function(array1, array2, flag) {
    difference <- sd(mapply('-', array1, array2, SIMPLIFY = TRUE), na.rm = TRUE)

    if (flag) {
        #calculate t_test for analysis genes
        if (difference != 0) {
            A <- t.test(array1, array2, var.equal = FALSE)[c('statistic',
                'p.value')]
            return(c(A$statistic, A$p.value))
        } else{
            return(c(NA, NA))
        }
    }
    else{
        #calculate Kolmogorov-Smirnov for cg analysis
        if (difference != 0 || is.na(difference)) {
            tryCatch({
                A <- ks.test(array1, array2)[c('p.value')]
            }, error = function(error_condition) { A <- list(p.value = NA) })
            return(c(A$p.value))
        }
    }
}
#[END]

setRowNames <- function(df) {
    rowNames <- c("medianDown","medianMedium","medianUP","bd_UPvsMID",
        "bd_UPvsDOWN","bd_MIDvsDOWN","pvalue_UPvsMID","pvalue_UPvsDOWN",
        "pvalue_MIDvsDOWN","pearson_correlation","pvalue_pearson_correlation",
        "fc_UPvsMID(gene)","fc_UPvsDOWN(gene)","fc_MIDvsDOWN(gene)",
        "pvalue_UPvsMID(gene)","pvalue_UPvsDOWN(gene)","pvalue_MIDvsDOWN(gene)"
    )

    for (i in 1:17) {
        row.names(df)[i] <- rowNames[i]
    }

    return(df)
}
