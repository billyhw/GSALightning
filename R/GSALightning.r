GSAfunc <- function(eset, fac, nperm, method = c('maxmean','mean','absmean')) {

    fac <- as.numeric(as.factor(fac))-1

    if (sum(fac==1) > sum(fac==0)) fac <- abs(fac - 1)

    numX <- sum(fac==1)
    numY <- sum(fac==0)

    permMat <- fac%*%t(rep(1,nperm))
    for (i in 1:nperm) permMat[,i] <- sample(permMat[,i])

    permMat <- Matrix(permMat)

    sumAll <- rowSums(eset)
    sumsqAll <- rowSums(eset^2)

    SumX <- eset %*% permMat
    SumXsq <- eset^2 %*% permMat

    SumY <- sumAll - SumX
    SumYsq <- sumsqAll - SumXsq

    MeanX <- SumX/numX
    MeanY <- SumY/numY

    resultsMat <- {(MeanX-MeanY)/sqrt(SumXsq-numX*MeanX^2+SumYsq-numY*MeanY^2)}*sqrt((length(fac)-2)/(1/numX+1/numY))
    if (method == 'mean') return(resultsMat)
    else if (method == 'absmean') return( abs(resultsMat) )
    else if (method == 'maxmean') {
        resultsMat1 <- pmax(as.matrix(resultsMat),0)
        resultsMat2 <- -pmin(as.matrix(resultsMat),0)
        return(ls=list(resultsMat1=resultsMat1,resultsMat2=resultsMat2))
    }
}

pvalFromPermMat <- function (obs, perms) {
    tempObs <- rep(obs, ncol(perms))
    dim(tempObs) <- dim(perms)
    rowSums(perms <= tempObs)
}

permTestLight <- function(eset, fac, nperm, tests = c('unpaired','paired','multi','wilcox'), method = c('mean','absmean'), npermBreaks = 2000, verbose = TRUE) {

    verbose <- isTRUE(verbose)

    tests <- match.arg(tests)
    method <- match.arg(method)

    mat <- Diagonal(nrow(eset),1)
    rownames(mat) <- colnames(mat) <- rownames(eset)

    # message("Note that permTestLight() simply runs GSALight() by treating each individual gene as a gene set.")

    results <- GSALight(eset, fac, mat, nperm, tests, method, restandardize = FALSE, npermBreaks=npermBreaks, verbose = verbose)
    results <- results[,-ncol(results)]

}

GSALight <- function (eset, fac, gs, nperm = NULL, tests = c('unpaired','paired','multi','wilcox'), method = c('maxmean','mean','absmean'),
                       minsize = 1, maxsize = Inf, restandardize = TRUE, npermBreaks = 2000,
                       rmGSGenes = c('stop', 'gene', 'gs'), verbose = TRUE) {

    restandardize <- isTRUE(restandardize)
    verbose <- isTRUE(verbose)

    tests <- match.arg(tests)
    method <- match.arg(method)
    rmGSGenes <- match.arg(rmGSGenes)

    if (tests == 'paired') {
        if (! class(fac) %in% c("numeric", "integer")) stop("For paired t-test, fac must be an integer vector.")
        if (any(fac == 0)) stop("For paired t-test, 0 is not allowed for subject pair index.")
        tab <- table(abs(fac))
        if (any(tab != 2) | sum(fac > 0) != sum(fac < 0)) stop("Some values in fac are not paired. For paired t-tests, fac must be a vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions.")
        sortedfac <- sort(unique(abs(fac)))
        if (!all(sortedfac == c(1:(length(fac)/2)))) stop("Some values in fac are skipped. For paired t-tests, the numbering of the subject pairs must be 1,2,3,... and -1,-2,-3,... with no skipped integers.")
    }

    if (tests == 'unpaired' | tests == 'wilcox') {
        if (! is.factor(fac)) fac <- as.factor(fac)
        if (length(unique(fac)) > 2 ) stop("More than two classes detected. For tests = 'unpaired' or 'wilcox' only two-sample tests are supported.")
    }

    if (tests == 'multi') {
        if (method != 'mean') stop("For the multi-class test, only the 'mean' method is available. Please set the method to 'mean' and restart.")
    }

    if (is.data.table(gs)) mat <- dataTable2Mat(gs)
    else if (is.list(gs)) {
        if (is.null(names(gs))) stop("Gene set names missing: each element in the gene set list must be named.")
        gs <- list2DataTable(gs)
        mat <- dataTable2Mat(gs)
    }
    else mat <- gs

    rm(gs)

    if (is.null(rownames(eset))) stop("Gene names missing: each row of the expression data must be named after the gene.")
    if (any(is.na(eset))) stop('The Expression data contain missing values.')
    if (any(apply(eset, 1, var) == 0)) stop('Some genes has 0 sample variance, please remove those genes prior to running GSALightning. Also consider removing genes with small sample variance.')

    if (! all(colnames(mat) %in% rownames(eset))) {
        if (rmGSGenes == 'gene') {
            if (verbose) message("Some genes within the gene sets are not contained in the expression data set.\n These genes are removed from the gene sets since rmGSGenes == 'gene'.")
            mat <- mat[,colnames(mat) %in% rownames(eset)]
        }
        else if (rmGSGenes == 'gs') {
            if (verbose) message("Some genes within the gene sets are not contained in the expression data set.\n Gene sets with missing genes are removed since rmGSGenes == 'gs'.")
            numGenes <- rowSums(mat)
            newNumGenes <- rowSums(mat[,colnames(mat) %in% rownames(eset)])
            mat <- mat[numGenes == newNumGenes,]
        }
        else stop("Some genes within the gene sets are not contained in the expression data set.\n Set rmGSGenes = 'gene' or 'gs' to remove respectively the missing genes or gene sets.")
    }

    setsize <- rowSums(mat)
    mat <- mat[setsize >= minsize & setsize <= maxsize,]

    mat <- mat[,colSums(mat) >= 1]
    eset <- eset[colnames(mat),]

    numGenes <- rowSums(mat)

    if (is.null(nperm)) {
        nperm <- nrow(mat)/0.05*2
        message("Number of permutations is not specified. Automatically set to ", nperm, ".")
    }

    if (verbose) message("After gene set size filtering, there are ", nrow(mat), " gene sets,\n containing a total of ", nrow(eset), " genes for analysis.")

    if (tests == 'wilcox') {
        message("Converting gene expression to ranks for Wilcoxon-Rank-Sum Test")
        for (i in 1:nrow(eset)) eset[i, ] <- rank(eset[i, ])
    }

    if (verbose) message("Obtaining observed gene set statistics.")

    if (tests == 'paired') obs <- rowPairedTtests(as.matrix(eset), fac, method)
    if (tests == 'multi') obs <- rowMultitests(as.matrix(eset), fac)
    if (tests == 'wilcox') obs <- rowWilcoxTest(as.matrix(eset), fac, method)
    if (tests == 'unpaired') obs <- rowtests(as.matrix(eset), fac, method)

    if (restandardize) {
        numCatGenes <- colSums(mat)
        totCatGenes <- sum(numCatGenes)
        if (method == 'maxmean') {
            meanobs1 <- sum(obs$results1*numCatGenes)/totCatGenes
            sdobs1 <- sqrt({sum(numCatGenes*{obs$results1^2}) - totCatGenes*meanobs1^2}/(totCatGenes - 1))
            meanobs2 <- sum(obs$results2*numCatGenes)/totCatGenes
            sdobs2 <- sqrt({sum(numCatGenes*{obs$results2^2}) - totCatGenes*meanobs2^2}/(totCatGenes - 1))
            obs1 <- as.vector(mat %*% obs$results1)
            obs2 <- as.vector(mat %*% obs$results2)
            obs1 <- obsOrig1 <- (obs1/numGenes - meanobs1)/sdobs1
            obs2 <- obsOrig2 <- (obs2/numGenes - meanobs2)/sdobs2
            obs <- pmax(obs1,obs2)
            obs[obs2 > obs1] <- -1*obs[obs2 > obs1]
            obsOrig <- pmax(obsOrig1,obsOrig2)
            obsOrig[obsOrig2 > obsOrig1] <- -1*obsOrig[obsOrig2 > obsOrig1]
        }
        else {
            meanobs <- sum(obs*numCatGenes)/totCatGenes
            sdobs <- sqrt({sum(numCatGenes*{obs^2}) - totCatGenes*meanobs^2}/(totCatGenes - 1))
            obs <- mat %*% obs
            obs <- as.vector(obs)
            obs <- obsOrig <- (obs/numGenes - meanobs)/sdobs*sqrt(numGenes)
        }
    }
    else {
        if (method == 'maxmean') {
            obs1 <- as.vector(mat %*% obs$results1)
            obs2 <- as.vector(mat %*% obs$results2)
            obsOrig1 <- obs1/numGenes
            obsOrig2 <- obs2/numGenes
            obs <- pmax(obs1,obs2)
            obs[obs2 > obs1] <- -1*obs[obs2 > obs1]
            obsOrig <- pmax(obsOrig1,obsOrig2)
            obsOrig[obsOrig2 > obsOrig1] <- -1*obsOrig[obsOrig2 > obsOrig1]
        }
        else {
            obs <- mat %*% obs
            obs <- as.vector(obs)
            obsOrig <- obs/numGenes
        }
    }

    if (nperm <= npermBreaks) {
        if (tests == 'paired') permMat <- GSApairedfunc(as.matrix(eset),fac,nperm,method)
        if (tests == 'multi') permMat <- GSAMultifunc(as.matrix(eset), fac, nperm)
        if (tests == 'wilcox') permMat <- GSAWilcoxFunc(as.matrix(eset), fac, nperm,method)
        if (tests == 'unpaired') permMat <- GSAfunc(as.matrix(eset),fac,nperm,method)
        if (restandardize) permMat <- getGeneSetPermStats(restandardize, method, permMat, mat, numCatGenes, totCatGenes, numGenes)
        else permMat <- getGeneSetPermStats(restandardize, method, permMat, mat)
        pvalSums <- pvalFromPermMat(obs, permMat)
        pval <- pvalSums/nperm
    }
    else {
        if (verbose) message("Running batch-mode permutation.")
        permBreaks <- ceiling(nperm/npermBreaks)
        pvalSums <- rep(0,nrow(mat))
        for (i in 1:permBreaks) {
            if (verbose) message("Running batch ", i, " of ", permBreaks, " batches.")
            if (tests == 'paired') permMat <- GSApairedfunc(as.matrix(eset),fac,ifelse(i!=permBreaks,npermBreaks,nperm-{i-1}*npermBreaks),method)
            if (tests == 'multi') permMat <- GSAMultifunc(as.matrix(eset),fac,ifelse(i!=permBreaks,npermBreaks,nperm-{i-1}*npermBreaks))
            if (tests == 'wilcox') permMat <- GSAWilcoxFunc(as.matrix(eset), fac, ifelse(i!=permBreaks, npermBreaks, nperm-{i-1}*npermBreaks), method)
            if (tests == 'unpaired') permMat <- GSAfunc(as.matrix(eset),fac,ifelse(i!=permBreaks,npermBreaks,nperm-{i-1}*npermBreaks),method)
            if (restandardize) permMat <- getGeneSetPermStats(restandardize, method, permMat, mat, numCatGenes, totCatGenes, numGenes)
            else permMat <- getGeneSetPermStats(restandardize, method, permMat, mat)
            pvalSums <- pvalSums + pvalFromPermMat(obs, permMat)
        }
        pval <- pvalSums/nperm
    }

    if (verbose) message("Permutation done. Evaluating P-values.")
    if (method == 'absmean' || tests == "multi") {
        pvals <- 1-pval
        qvals <- p.adjust(pvals,method='BH')
        results <- cbind(pvals,qvals,obsOrig,rowSums(mat))
        colnames(results) <- c('p-value','q-value','statistics','# genes')
    }
    else {
        if (tests == 'unpaired' || tests == 'wilcox') {
            lvls <- levels(as.factor(fac))
            fac <- as.numeric(as.factor(fac))-1
            if (sum(fac==1) > sum(fac==0)) lvls <- rev(lvls)

            pvals <- cbind(pval,1-pval)
            rownames(pvals) <- rownames(mat)

            qvals <- cbind(p.adjust(pvals[,1],'BH'),p.adjust(pvals[,2],'BH'))

            results <- cbind(pvals,qvals,obsOrig,rowSums(mat))
            colnames(results) <- c(paste('p-value:up-regulated in', lvls[1]), paste('p-value:up-regulated in',lvls[2]),
                                   paste('q-value:up-regulated in', lvls[1]), paste('q-value:up-regulated in',lvls[2]),
                                   paste('statistics: up-regulated in',lvls[2]),'# genes')
        }
        else {
            pvals <- cbind(1-pval,pval)
            qvals <- cbind(p.adjust(pvals[,1],'BH'),p.adjust(pvals[,2],'BH'))

            results <- cbind(pvals,qvals,obsOrig,rowSums(mat))
            colnames(results) <- c(paste('p-value:up-regulated in positives'), paste('p-value:up-regulated in negatives'),
                                   paste('q-value:up-regulated in positives'), paste('q-value:up-regulated in negatives'),
                                   'statistics (up-regulated in positives)','# genes')
        }
    }

    as.data.frame(results)

}

getGeneSetPermStats <- function(restandardize, method, permMat, mat, numCatGenes=NULL, totCatGenes=NULL, numGenes=NULL) {
    if (restandardize) {
        if (method == 'maxmean') {
            meanStar1 <- colSums(permMat$resultsMat1*numCatGenes)/totCatGenes
            sdStar1 <- sqrt((colSums({permMat$resultsMat1^2}*numCatGenes) - totCatGenes*meanStar1^2)/(totCatGenes - 1))
            permMat1 <- as.matrix(mat%*%permMat$resultsMat1)
            permMat1 <- t((t(permMat1/numGenes) - meanStar1)/sdStar1)

            meanStar2 <- colSums(permMat$resultsMat2*numCatGenes)/totCatGenes
            sdStar2 <- sqrt((colSums({permMat$resultsMat2^2}*numCatGenes) - totCatGenes*meanStar2^2)/(totCatGenes - 1))
            permMat2 <- as.matrix(mat%*%permMat$resultsMat2)
            permMat2 <- t((t(permMat2/numGenes) - meanStar2)/sdStar2)

            permMat <- pmax(permMat1, permMat2)
            permMat[permMat2 > permMat1] <- -1 * permMat[permMat2 > permMat1]
        }
        else {
            meanStar <- colSums(permMat*numCatGenes)/totCatGenes
            sdStar <- sqrt((colSums({permMat^2}*numCatGenes) - totCatGenes*meanStar^2)/(totCatGenes - 1))
            permMat <- as.matrix(mat%*%permMat)
            permMat <- t((t(permMat/numGenes) - meanStar)/sdStar)*sqrt(numGenes)
        }
    }
    else {
        if (method == 'maxmean') {
            permMat1 <- as.matrix(mat%*%permMat$resultsMat1)
            permMat2 <- as.matrix(mat%*%permMat$resultsMat2)
            permMat <- pmax(permMat1, permMat2)
            permMat[permMat2 > permMat1] <- -1 * permMat[permMat2 > permMat1]
        }
        else permMat <- as.matrix(mat%*%permMat)
    }
    permMat
}

rowtests <- function(eset,fac,method=c('maxmean','mean','absmean')) {
    fac <- as.numeric(as.factor(fac))-1
    if (sum(fac==1) > sum(fac==0)) fac <- abs(fac - 1)
    numx <- sum(fac==1)
    numy <- sum(fac==0)
    x1 <- eset[,fac==1]
    x2 <- eset[,fac==0]
    mean1 <- rowMeans(x1)
    mean2 <- rowMeans(x2)
    var1 <- rowSums(x1^2) - numx*mean1^2
    var2 <- rowSums(x2^2) - numy*mean2^2
    results <- {mean1-mean2}/sqrt({{var1 + var2}/{length(fac)-2}}*{1/numx+1/numy})
    if (method == 'mean') return( results )
    else if (method == 'absmean') return( abs(results) )
    else if (method == 'maxmean') {
        results1 <- pmax(results,0)
        results2 <- -pmin(results,0)
        return(ls=list(results1=results1,results2=results2))
    }
}

rowPairedTtests <- function(eset, fac, method = c('maxmean','mean','absmean')) {

    numSubject <- length(fac)/2
    mat <- Matrix(0,ncol(eset),numSubject)
    for (i in 1:numSubject) mat[which(abs(fac) == i),i] <- sign(fac[which(abs(fac)==i)])

    d <- eset%*%mat
    meand <- rowMeans(d)
    vard <- rowSums((d - meand)^2)/(numSubject-1)
    results <- meand/(sqrt(vard)/sqrt(numSubject))

    if (method == 'mean') return(results)
    else if (method == 'absmean') return( abs(results) )
    else if (method == 'maxmean') {
        results1 <- pmax(results,0)
        results2 <- -pmin(results,0)
        return(ls=list(results1=results1,results2=results2))
    }
}

GSApairedfunc <- function(eset, fac, nperm, method = c('maxmean','mean','absmean')) {

    numSubject <- length(fac)/2
    mat <- Matrix(0,ncol(eset),numSubject)
    for (i in 1:numSubject) mat[which(abs(fac) == i),i] <- sign(fac[which(abs(fac)==i)])

    d <- eset%*%mat

    permMat <- matrix(0,numSubject,nperm)

    for (i in 1:nperm) permMat[,i] <- rbinom(numSubject,1,prob=0.5)

    permMat <- Matrix(permMat)

    sumd <- rowSums(d)
    sumpermd <- 2*d%*%permMat
    meanx <- (sumd-sumpermd)/numSubject
    sumdsq <- rowSums(d^2)
    sdx <- sqrt((sumdsq-numSubject*meanx^2)/(numSubject-1))
    resultsMat <- meanx/(sdx/sqrt(numSubject))

    if (method == 'mean') return(resultsMat)
    else if (method == 'absmean') return( abs(resultsMat) )
    else if (method == 'maxmean') {
        resultsMat1 <- pmax(as.matrix(resultsMat),0)
        resultsMat2 <- -pmin(as.matrix(resultsMat),0)
        return(ls=list(resultsMat1=resultsMat1,resultsMat2=resultsMat2))
        }
}

list2DataTable <- function(geneList) {
    gsLength <- sapply(geneList,length)
    data.table(geneSet = rep(names(geneList),times=gsLength), gene = unlist(geneList))
}

dataTable2Mat <- function(gsTable) {
    tmp <- factor(gsTable$gene)
    tmpGenes <- levels(tmp)
    geneFactor <- as.numeric(tmp)
    origOrder <- unique(gsTable$geneSet)
    tmp <- factor(gsTable$geneSet)
    tmpGS <- levels(tmp)
    gs <- as.numeric(tmp)

    mat <- sparseMatrix(gs,geneFactor,x=1)

    if (any(mat > 1)) {
        warning("There appears to be duplicated genes in some gene sets. The duplicated genes are removed.")
        mat[mat > 1] <- 1
    }
    rownames(mat) <- tmpGS
    colnames(mat) <- tmpGenes

    mat[origOrder,]
}

rowMultitests <- function(eset,fac) {

    fac <- as.numeric(as.factor(fac))-1

    #if (sum(fac==1) > sum(fac==0)) fac <- abs(fac - 1)

    numPerClass <- table(fac)

    totalMean <- rowMeans(eset)

    classMean <- withinClassSS <- matrix(0, nrow(eset), length(numPerClass))

    for (i in 0:(length(numPerClass)-1)) {
        classMean[, i+1] <- rowMeans(eset[, fac == i])
        withinClassSS[, i+1] <- rowSums((eset[, fac == i] - classMean[, i+1])^2)
    }

    numerMult <- length(fac)/prod(numPerClass)
    denomMult <- sum(1/numPerClass)/sum(numPerClass-1)

    numer <- sqrt(numerMult*colSums(t(classMean - totalMean)^2*as.vector(numPerClass)))
    denom <- sqrt(denomMult*rowSums(withinClassSS))

    numer/denom

}

GSAMultifunc <- function(eset, fac, nperm) {

    # turn factor to numeric
    fac <- as.numeric(as.factor(fac))-1

    # get number of subject per class
    numPerClass <- table(fac)
    #numX <- sum(fac==1)
    #numY <- sum(fac==0)

    permMat <- fac%*%t(rep(1,nperm))
    for (i in 1:nperm) permMat[,i] <- sample(permMat[,i])

    permMat <- Matrix(permMat)

    meanAll <- rowMeans(eset)

    meanMatList <- vector("list", length(numPerClass))
    meanMatList <- lapply(meanMatList, function (x) matrix(0, nrow(eset), nperm))

    # SSWithinClass <- matrix(0, nrow(eset), nperm)
    for (i in 0:(length(numPerClass)-1)) {

        tmpPermMat <- Matrix(0, length(fac), nperm)
        # tmpPermMat[permMat != i] <- 0
        tmpPermMat[permMat == i] <- 1

        meanMatList[[i+1]] <- (eset %*% tmpPermMat)/numPerClass[i+1]
        if (i == 0) {
            SSBetweenClass <- numPerClass[1]*(meanMatList[[1]] - meanAll)^2
            SSWithinClass <- eset^2 %*% tmpPermMat - numPerClass[i+1]*meanMatList[[i+1]]^2
        }
        else {
            SSBetweenClass <- SSBetweenClass + numPerClass[i+1]*(meanMatList[[i+1]] - meanAll)^2
            SSWithinClass <- SSWithinClass + eset^2 %*% tmpPermMat - numPerClass[i+1]*meanMatList[[i+1]]^2
        }
    }

    numerMult <- length(fac)/prod(numPerClass)
    denomMult <- sum(1/numPerClass)/sum(numPerClass-1)

    numer <- sqrt(numerMult*SSBetweenClass)
    denom <- sqrt(denomMult*SSWithinClass)

    numer/denom

}

rowWilcoxTest <- function(esetRank, fac, method = c('maxmean','mean','absmean')) {
    fac <- as.numeric(as.factor(fac))-1
    if (sum(fac==1) > sum(fac==0)) fac <- abs(fac - 1)
    numX <- sum(fac==1)
    numY <- sum(fac==0)
    results <- esetRank%*%fac - numX*(numX + numY + 1)/2
    results <- results/sqrt(numX*numY*(numX + numY + 1)/12)
    if (method == 'mean') return( results )
    else if (method == 'absmean') return( abs(results) )
    else if (method == 'maxmean') {
        results1 <- pmax(results,0)
        results2 <- -pmin(results,0)
        return(ls=list(results1=results1,results2=results2))
    }
}

GSAWilcoxFunc <- function(esetRank, fac, nperm, method = c('maxmean','mean','absmean')) {

    fac <- as.numeric(as.factor(fac))-1

    if (sum(fac==1) > sum(fac==0)) fac <- abs(fac - 1)

    numX <- sum(fac==1)
    numY <- sum(fac==0)

    permMat <- fac%*%t(rep(1,nperm))
    for (i in 1:nperm) permMat[,i] <- sample(permMat[,i])

    permMat <- Matrix(permMat)

    resultsMat <- esetRank %*% permMat - numX*(numX + numY + 1)/2

    resultsMat <- resultsMat/sqrt(numX*numY*(numX + numY + 1)/12)

    if (method == 'mean') return(resultsMat)
    else if (method == 'absmean') return( abs(resultsMat) )
    else if (method == 'maxmean') {
        resultsMat1 <- pmax(as.matrix(resultsMat),0)
        resultsMat2 <- -pmin(as.matrix(resultsMat),0)
        return(ls=list(resultsMat1=resultsMat1,resultsMat2=resultsMat2))
    }
}
