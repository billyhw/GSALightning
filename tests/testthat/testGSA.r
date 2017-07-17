# to install the testing branch of GSA-Lightning
# library(githubinstall)
# gh_install_packages("GSALighting", ref = "testing")

# to run, type:
# test_dir("tianLab/GSALightning/tests/testthat")

library(GSA)
library(GSALightning)

set.seed(100)
x <- matrix(rnorm(1000*20),ncol=20)
rownames(x) <- paste("g",1:1000,sep="")
dd <- sample(1:1000,size=100)

genenames=paste("g",1:1000,sep="")
y <- factor(c(rep(0,10),rep(1,10)))
yy <- factor(c(rep(1,10),rep(2,10)))

#create some random gene sets
genesets=vector("list",50)
for(i in 1:50){
 genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
}
names(genesets)=paste("set",as.character(1:50),sep="")

test_that("Unpaired Test: Concordance of Gene set Statistics and p-values with GSA", {
  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'unpaired',
    method = 'mean', restandardize = TRUE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "mean",
    resp.type = "Two class unpaired", s0 = 0, s0.perc=-1, restand = TRUE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,2], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'unpaired',
    method = 'mean', restandardize = FALSE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "mean",
    resp.type = "Two class unpaired", s0 = 0, s0.perc=-1, restand = FALSE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,2], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'unpaired',
    method = 'absmean', restandardize = TRUE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "absmean",
    resp.type = "Two class unpaired", s0 = 0, s0.perc=-1, restand = TRUE)

  #expect_equivalent(GSALightResult[,3], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'unpaired',
    method = 'absmean', restandardize = FALSE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "absmean",
    resp.type = "Two class unpaired", s0 = 0, s0.perc=-1, restand = FALSE)

  expect_equivalent(GSALightResult[,3], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'unpaired',
    method = 'maxmean', restandardize = TRUE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "maxmean",
    resp.type = "Two class unpaired", s0 = 0, s0.perc=-1, restand = TRUE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,2], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'unpaired',
    method = 'maxmean', restandardize = FALSE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "maxmean",
    resp.type = "Two class unpaired", s0 = 0, s0.perc=-1, restand = FALSE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,2], GSAResult$pvalues.hi), 0.90)
})

# classes for paired-data
y <- yy <- c(1:10, -(1:10))

test_that("Paired-Test: Concordance of Gene set Statistics and p-values with GSA", {
  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'paired',
    method = 'mean', restandardize = TRUE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "mean",
    resp.type = "Two class paired", s0 = 0, s0.perc=-1, restand = TRUE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'paired',
    method = 'mean', restandardize = FALSE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "mean",
    resp.type = "Two class paired", s0 = 0, s0.perc=-1, restand = FALSE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'paired',
    method = 'absmean', restandardize = TRUE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "absmean",
    resp.type = "Two class paired", s0 = 0, s0.perc=-1, restand = TRUE)

  #expect_equivalent(GSALightResult[,3], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'paired',
    method = 'absmean', restandardize = FALSE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "absmean",
    resp.type = "Two class paired", s0 = 0, s0.perc=-1, restand = FALSE)

  expect_equivalent(GSALightResult[,3], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'paired',
    method = 'maxmean', restandardize = TRUE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "maxmean",
    resp.type = "Two class paired", s0 = 0, s0.perc=-1, restand = TRUE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)

  GSALightResult <- GSALight(x, y, genesets, nperm = 200, tests = 'paired',
    method = 'maxmean', restandardize = FALSE)
  GSAResult <- GSA(x, yy, genesets, genenames, nperms=200, method = "maxmean",
    resp.type = "Two class paired", s0 = 0, s0.perc=-1, restand = FALSE)

  expect_equivalent(GSALightResult[,5], GSAResult$GSA.scores)
  expect_gt(cor(GSALightResult[,1], GSAResult$pvalues.hi), 0.90)
})

