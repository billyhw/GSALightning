# to install the testing branch of GSA-Lightning
# library(githubinstall)
# gh_install_packages("GSALighting", ref = "testing")

# to check coverage use package_coverage() in package "covr"

set.seed(100)
x <- matrix(rnorm(1000*20),ncol=20)
rownames(x) <- paste("g",1:1000,sep="")
dd <- sample(1:1000,size=100)

genenames=paste("g",1:1000,sep="")

y <- factor(c(rep(0,10),rep(1,10)))
yy <- c(1:10, -(1:10))
yyy <- factor(c(rep(0,7),rep(1,7),rep(2,6)))

#create some random gene sets
genesets=vector("list",50)
for(i in 1:50){
 genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
}
names(genesets)=paste("set",as.character(1:50),sep="")

test_that("Unpaired test only accept 2-classes", {
  expect_error(GSALight(x, yyy, genesets),
    "More than two classes detected. GSALightning only supports two-sample t-tests.")
})

test_that("Ensure paired-test class-labels correctness", {

  expect_error(GSALight(x, as.factor(y), genesets, tests = "paired"),
    "For paired t-test, fac must be an integer vector.")

  badYY <- c(0:9, -(0:9))
  expect_error(GSALight(x, badYY, genesets, tests = "paired"),
    "For paired t-test, 0 is not allowed for subject pair index.")

  badYY <- c(1:10, -(1:9), 100)
  expect_error(GSALight(x, badYY, genesets, tests = "paired"),
    "Some values in fac are not paired. For paired t-tests, fac must be a vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions.")

  badYY <- c(-100, 1:9, -(1:9), 100)
  expect_error(GSALight(x, badYY, genesets, tests = "paired"),
    "Some values in fac are skipped. For paired t-tests, the numbering of the subject pairs must be 1,2,3,... and -1,-2,-3,... with no skipped integers.")

})

test_that("Multiclass test must use 'mean' method", {
  expect_error(GSALight(x, yyy, genesets, tests = "multi", method = "absmean"),
    "For the multi-class test, only the 'mean' method is available. Please set the method to 'mean' and restart.")
})

test_that("Gene set must be named", {
  badGenesets = genesets
  names(badGenesets) = NULL
  expect_error(GSALight(x, y, badGenesets),
    "Gene set names missing: each element in the gene set list must be named.")
})

test_that("Check automatic number of permutations.", {
  nperm <- length(genesets)/0.05*2
  expect_message(GSALight(x, y, genesets, verbose = FALSE),
    paste("Number of permutations is not specified. Automatically set to ", nperm, ".\n", sep = ""))
  })

test_that("Gene names as rows in expression matrix", {
  badX <- x
  rownames(badX) <- NULL
  expect_error(GSALight(badX, y, genesets),
    "Gene names missing: each row of the expression data must be named after the gene.")
  })

test_that("Expression matrix cannot contain missing values", {
  badX <- x
  badX[1,1] <- NA
  expect_error(GSALight(badX, y, genesets),
    "The Expression data contain missing values.")
  })

test_that("Genes with 0 sample variance", {
  badX <- x
  badX[1,] <- 10
  expect_error(GSALight(badX, y, genesets),
    "Some genes has 0 sample variance, please remove those genes prior to running GSALightning. Also consider removing genes with small sample variance.")
  })

    # if (! all(colnames(mat) %in% rownames(eset))) {
    #     if (rmGSGenes == 'gene') {
    #         if (verbose) message("Some genes within the gene sets are not contained in the expression data set.\n These genes are removed from the gene sets since rmGSGenes == 'gene'.")
    #         mat <- mat[,colnames(mat) %in% rownames(eset)]
    #     }
    #     else if (rmGSGenes == 'gs') {
    #         if (verbose) message("Some genes within the gene sets are not contained in the expression data set.\n Gene sets with missing genes are removed since rmGSGenes == 'gs'.")
    #         numGenes <- rowSums(mat)
    #         newNumGenes <- rowSums(mat[,colnames(mat) %in% rownames(eset)])
    #         mat <- mat[numGenes == newNumGenes,]
    #     }
    #     else stop("Some genes within the gene sets are not contained in the expression data set.\n Set rmGSGenes = 'gene' or 'gs' to remove respectively the missing genes or gene sets.")
    # }

test_that("Error of gene sets contain genes not in expression set", {
  badGenesets = genesets
  badGenesets[[1]] = c(badGenesets[[1]], "aBadGene")
  expect_error(GSALight(x, y, badGenesets),
    "Some genes within the gene sets are not contained in the expression data set.\n Set rmGSGenes = 'gene' or 'gs' to remove respectively the missing genes or gene sets.")
})

test_that("Remove gene sets containing genes not in expression set (rmGSGenes = 'gs')", {
  badGenesets = genesets
  badGenesets[[1]] = c(badGenesets[[1]], "aBadGene")
  expect_message(GSALight(x, y, badGenesets, rmGSGenes = "gs"),
    "Some genes within the gene sets are not contained in the expression data set.\n Gene sets with missing genes are removed since rmGSGenes == 'gs'.")
  result = GSALight(x, y, badGenesets, rmGSGenes = "gs")
  expect_equal(rownames(result), names(genesets)[-1])
})

test_that("Remove genes not in expression set from gene sets (rmGSGenes = 'gene')", {
  badX = x
  badX = badX[-1,]
  expect_message(GSALight(x, y, badGenesets, rmGSGenes = "gene"),
    "Some genes within the gene sets are not contained in the expression data set.\n These genes are removed from the gene sets since rmGSGenes == 'gene'.")
  result = GSALight(badX, y, genesets, rmGSGenes = "gene")
  ind = sapply(genesets, function(x) "g1" %in% x)
  expect_equal(result[ind,6], rep(29, sum(ind)))
  expect_equal(result[!ind,6], rep(30, sum(!ind)))
})

