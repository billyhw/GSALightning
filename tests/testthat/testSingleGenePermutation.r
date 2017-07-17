# to install the testing branch of GSA-Lightning
# library(githubinstall)
# gh_install_packages("GSALighting", ref = "testing")

library(samr)
library(GSALightning)

set.seed(100)
x <- matrix(rnorm(1000*20),ncol=20)
rownames(x) <- paste("g",1:1000,sep="")
dd <- sample(1:1000,size=100)

genenames=paste("g",1:1000,sep="")
y <- factor(c(rep(0,10),rep(1,10)))

data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=genenames, logged2=FALSE)

# function to match p-values with SAM
matchPvals = function(GSAResult) {
  GSALightPvals <- rep(NA, nrow(GSAResult))
  for (i in 1:nrow(GSAResult)) {
    if (GSAResult[i,5] >= 0) GSALightPvals[i] <- GSAResult[i,2]
    else GSALightPvals[i] <- GSAResult[i,1]
  }
  GSALightPvals
}

test_that("Unpaired Test: Concordance of Statistics and p-values with SAM", {

  samr.obj<-samr(data, resp.type="Two class unpaired", nperms=1000, testStatistic=c("standard"), s0.perc=-1)
  GSAResult <- permTestLight(x, y, nperm = 1000, tests = 'unpaired', method = 'mean')

  expect_equivalent(GSAResult[,5], samr.obj$numer/(samr.obj$sd - samr.obj$s0))

  samrPVals <- samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
  expect_gt(cor(samrPVals, matchPvals(GSAResult)), 0.98)

})

test_that("Wilcoxon Rank Sum Test: Concordance of Statistics and p-values with SAM", {

  samr.obj<-samr(data, resp.type="Two class unpaired", nperms=200, testStatistic=c("wilcox"), s0.perc=-1)
  GSAResult <- permTestLight(x, y, nperm = 200, tests = 'wilcox', method = 'mean')

  expect_equivalent(GSAResult[,5], samr.obj$numer/(samr.obj$sd - samr.obj$s0))

  samrPVals <- samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
  expect_gt(cor(samrPVals, matchPvals(GSAResult)), 0.95)

})

# classes for paired-data
y <- c(1:10, -(1:10))

data=list(x=x,y=y, geneid=as.character(1:nrow(x)),
genenames=genenames, logged2=FALSE)

# function to match p-values with SAM
matchPvals = function(GSAResult) {
  GSALightPvals <- rep(NA, nrow(GSAResult))
  for (i in 1:nrow(GSAResult)) {
    if (GSAResult[i,5] >= 0) GSALightPvals[i] <- GSAResult[i,1]
    else GSALightPvals[i] <- GSAResult[i,2]
  }
  GSALightPvals
}

test_that("Paired-Test: Concordance of Gene set Statistics and p-values with GSA", {

  samr.obj<-samr(data, resp.type="Two class paired", nperms=1000, testStatistic=c("standard"), s0.perc=-1)
  GSAResult <- permTestLight(x, y, nperm = 1000, tests = 'paired', method = 'mean')

  expect_equivalent(GSAResult[,5], samr.obj$numer/(samr.obj$sd - samr.obj$s0))

  samrPVals <- samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
  expect_gt(cor(samrPVals, matchPvals(GSAResult)), 0.98)

})


