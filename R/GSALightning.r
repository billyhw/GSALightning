GSAfunc <- function(eset, fac, nperm, method = c('mean','absmean')) {
    
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
    
}

pvalFromPermMat <- function (obs, perms) {
    tempObs <- rep(obs, ncol(perms))
    dim(tempObs) <- dim(perms)
    rowSums(perms <= tempObs)
}

permTestLight <- function(eset, fac, nperm, tests = c('unpaired','paired'), method = c('mean','absmean'), npermBreaks = 2000, verbose = TRUE) {

    mat <- Diagonal(nrow(eset),1)
    rownames(mat) <- colnames(mat) <- rownames(eset)
    results <- GSALight(eset, fac, mat, nperm, tests, method, restandardize = F, npermBreaks=npermBreaks, verbose = verbose)
    results <- results[,-ncol(results)]
    
}

GSALight <- function (eset, fac, gs, nperm = NULL, tests = c('unpaired','paired'), method = c('mean','absmean'),
                      minsize = 1, maxsize = Inf, restandardize = T, npermBreaks = 2000,
                      rmGSGenes = c('stop', 'gene', 'gs'), verbose = TRUE) {

    if (length(tests) > 1) tests <- 'unpaired'
    if (length(method) > 1) method <- 'mean'
    if (length(rmGSGenes) > 1) rmGSGenes <- 'stop'
    if (! tests %in% c('unpaired','paired')) stop("tests must be one of 'unpaired' or 'paired'")
    if (! method %in% c('mean','absmean')) stop("method must be one of 'mean' or 'absmean'")
    if (!rmGSGenes %in% c('stop','gene','gs')) stop("rmGSGenes must be one of 'stop', 'gene', or 'gs'.")

    if (tests == 'paired') {
        if (! is.integer(fac)) stop(" For paired t-test, fac must be an integer vector.")
        if (any(fac == 0)) stop(" For paired t-test, 0 is not allowed for subject pair index.")
        tab <- table(abs(fac))
        if (any(tab != 2) | sum(fac > 0) != sum(fac < 0)) stop(" Some values in fac are not paired. For paired t-tests, fac must be a vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions. ")
        sortedfac <- sort(unique(abs(fac)))
        if (!all(sortedfac == c(1:(length(fac)/2)))) stop(" Some values in fac are skipped. For paired t-tests, the numbering of the subject pairs must be 1,2,3,..., with no skipped integers. ")
    }

    if (tests == 'unpaired') {
        if (! is.factor(fac)) fac <- as.factor(fac)
        if (length(unique(fac)) > 2 ) stop("more than two classes detected. GSALightning only supports two-sample t-tests.")
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
    if (any(is.na(eset))) stop('The Expression data contain missing values. \n')
    if (any(apply(eset, 1, var) == 0)) stop('Some genes has 0 sample variance, please remove those genes prior to running GSALightning. Also consider removing genes with small sample variance.')

    if (! all(colnames(mat) %in% rownames(eset))) {
        if (rmGSGenes == 'gene') {
            if (verbose) cat("Some genes within the gene sets are not contained in the expression data set. \n These genes are removed from the gene sets since rmGSGenes == 'gene'.\n")
            mat <- mat[,colnames(mat) %in% rownames(eset)]
        }
        else if (rmGSGenes == 'gs') {
            if (verbose) cat("Some genes within the gene sets are not contained in the expression data set. \n Gene sets with missing genes are removed since rmGSGenes == 'gs'.\n")
            numGenes <- rowSums(mat)
            newNumGenes <- rowSums(mat[,colnames(mat) %in% rownames(eset)])
            mat <- mat[numGenes == newNumGenes,]
        }
        else stop("Some genes within the gene sets are not contained in the expression data set. \n Set rmGSGenes = 'gene' or 'gs' to remove respectively the missing genes or gene sets.")
    }
    
    setsize <- rowSums(mat)
    mat <- mat[setsize >= minsize & setsize <= maxsize,]

    mat <- mat[,colSums(mat) >= 1]
    eset <- eset[colnames(mat),]

    numGenes <- rowSums(mat)
    
    if (is.null(nperm)) {
        nperm <- nrow(mat)/0.05*2
        cat("Number of permutations is not specified. Automatically set to", nperm, ". \n")
    }

    if (verbose) cat("After gene set size filtering, there are", nrow(mat), "gene sets, \n containing a total of", nrow(eset), "genes for analysis. \n")

    if (verbose) cat("Obtaining observed gene set statistics. \n")

    if (tests == 'paired') obs <- rowPairedTtests(as.matrix(eset), fac, method)
    else obs <- rowtests(as.matrix(eset), fac, method)

    if (restandardize) {
        numCatGenes <- colSums(mat)
        totCatGenes <- sum(numCatGenes)
        meanobs <- sum(obs*numCatGenes)/totCatGenes
        #sdobs <- sd(obs)
        sdobs <- sqrt({sum(numCatGenes*{obs^2}) - totCatGenes*meanobs^2}/(totCatGenes - 1))
    }
    
    obs <- mat %*% obs
    obs <- as.vector(obs)
    obsOrig <- obs/numGenes

    if (restandardize) obs <- (obs/numGenes - meanobs)/sdobs
    
    if (nperm <= npermBreaks) {
        if (verbose) cat("Running permutation. \n")
        if (tests == 'paired') permMat <- GSApairedfunc(as.matrix(eset),fac,nperm,method)
        else permMat <- GSAfunc(as.matrix(eset),fac,nperm,method)
        if (restandardize) {
            #meanStar <- colMeans(permMat)
            #sdStar <- sqrt((colSums(permMat^2) - nrow(eset)*meanStar^2)/(nrow(eset)-1))
            meanStar <- colSums(permMat*numCatGenes)/totCatGenes
            sdStar <- sqrt((colSums({permMat^2}*numCatGenes) - totCatGenes*meanStar^2)/(totCatGenes - 1))
            permMat <- as.matrix(mat%*%permMat)
            permMat <- t((t(permMat/numGenes) - meanStar)/sdStar)
        }
        else permMat <- as.matrix(mat%*%permMat)
        pvalSums <- pvalFromPermMat(obs, permMat)
        pval <- pvalSums/nperm
    }
    else {
        if (verbose) cat("Running batch-mode permutation. \n")
        permBreaks <- ceiling(nperm/npermBreaks)
        pvalSums <- rep(0,nrow(mat))
        for (i in 1:permBreaks) {
            if (verbose) cat("Running batch", i, "of", permBreaks, "batches. \n")
            if (tests == 'paired') permMat <- GSApairedfunc(as.matrix(eset),fac,ifelse(i!=permBreaks,npermBreaks,nperm-{i-1}*npermBreaks),method)
            else permMat <- GSAfunc(as.matrix(eset),fac,ifelse(i!=permBreaks,npermBreaks,nperm-{i-1}*npermBreaks),method)
            if (restandardize) {
                #meanStar <- colMeans(permMat)
                #sdStar <- sqrt((colSums(permMat^2) - nrow(eset)*meanStar^2)/(nrow(eset)-1))
                meanStar <- colSums(permMat*numCatGenes)/totCatGenes
                sdStar <- sqrt((colSums({permMat^2}*numCatGenes) - totCatGenes*meanStar^2)/(totCatGenes - 1))
                permMat <- as.matrix(mat%*%permMat)
                permMat <- t((t(permMat/numGenes) - meanStar)/sdStar)
            }
            else permMat <- as.matrix(mat%*%permMat)
            pvalSums <- pvalSums + pvalFromPermMat(obs, permMat)
        }
        pval <- pvalSums/nperm
    }

    if (verbose) cat("Permutation done. Evaluating P-values. \n")
    if (method == 'absmean') {
        pvals <- 1-pval
        qvals <- p.adjust(pvals,method='BH')
        results <- cbind(pvals,qvals,obsOrig,rowSums(mat))
        colnames(results) <- c('p-value','q-value','statistics','# genes')
    }
    else {
        if (tests == 'unpaired') {
            lvls <- levels(as.factor(fac))
            fac <- as.numeric(as.factor(fac))-1
            if (sum(fac==1) > sum(fac==0)) lvls <- rev(lvls)
        
            pvals <- cbind(1-pval,pval)
            colnames(pvals) <- c(paste('up-regulated in', lvls[2]), paste('up-regulated in',lvls[1]))
            rownames(pvals) <- rownames(mat)
        
            if (sum(fac==1) > sum(fac==0)) pvals <- pvals[,2:1]
            qvals <- cbind(p.adjust(pvals[,1],'BH'),p.adjust(pvals[,2],'BH'))

            results <- cbind(pvals,qvals,obsOrig,rowSums(mat))
            colnames(results) <- c(paste('p-value:up-regulated in', lvls[1]), paste('p-value:up-regulated in',lvls[2]),
                                   paste('q-value:up-regulated in', lvls[1]), paste('q-value:up-regulated in',lvls[2]),
                                   'statistics','# genes')
        }
        else {
            pvals <- cbind(1-pval,pval)
            qvals <- cbind(p.adjust(pvals[,1],'BH'),p.adjust(pvals[,2],'BH'))

            results <- cbind(pvals,qvals,obsOrig,rowSums(mat))
            colnames(results) <- c(paste('p-value:up-regulated in positives'), paste('p-value:up-regulated in negatives'),
                                   paste('q-value:up-regulated in positives'), paste('q-value:up-regulated in negatives'),
                                   'statistics','# genes')
        }
    }
    results
    
}

rowtests <- function(eset,fac,method=c('mean','absmean')) {
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
}

rowPairedTtests <- function(eset, fac, method = c('mean','absmean')) {

    numSubject <- length(fac)/2
    mat <- Matrix(0,ncol(eset),numSubject)
    for (i in 1:numSubject) mat[which(abs(fac) == i),i] <- sign(fac[which(abs(fac)==i)])

    d <- eset%*%mat
    meand <- rowMeans(d)
    vard <- rowSums((d - meand)^2)/(numSubject-1)
    results <- meand/(sqrt(vard)/sqrt(numSubject))
    
    if (method == 'mean') return(results)
    else if (method == 'absmean') return( abs(results) )

}

GSApairedfunc <- function(eset, fac, nperm, method = c('mean','absmean')) {

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

}

list2DataTable <- function(geneList) {
    gsLength <- sapply(geneList,length)
    data.table(geneSet = rep(names(geneList),times=gsLength), gene = unlist(geneList))
}

dataTable2Mat <- function(gsTable) {
    tmp <- factor(gsTable$gene)
    tmpGenes <- levels(tmp)
    geneFactor <- as.numeric(tmp)
    tmp <- factor(gsTable$geneSet)
    tmpGS <- levels(tmp)
    gs <- as.numeric(tmp)

    mat <- sparseMatrix(gs,geneFactor,x=1)

    rownames(mat) <- tmpGS
    colnames(mat) <- tmpGenes

    mat
}

