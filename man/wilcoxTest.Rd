\name{wilcoxTest}
\alias{wilcoxTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single-Gene Mann Whitney Wilcoxon Test
}
\description{
A two-sample single-gene differential expression analysis using the Mann Whitney Wilcoxon Test.
}
\usage{
wilcoxTest(eset, fac, tests = c("unpaired", "paired"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
The expression matrix. Each row is a gene, and each column is a subject/sample. The gene names must be presented as the row names.
}
  \item{fac}{
Subject labels, for unpaired T-test, either a factor or something that can be coerced into a factor (e.g. 0 and 1, Experiment and Control). For paired T-test,  fac must be an integer vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions.}
  \item{tests}{
The tests to performed. Can be either the default "unpaired" for unpaired T-tests or "paired" for paired T-tests.
}
}
\details{
This function performs Mann Whitney Wilcoxon test (a.k.a Mann Whitney U test and Wilcoxon Rank Sum test) for all genes in eset. This function is built on the wilcox.test() function in the "stats" package, but is structured to align with the usage of permTestLight. We included this function in case users may want to compare the permutation test results with the results of a standard, non-permutation-based nonparametric test.
}
\value{
A data frame with the p-values, the q-values (via Benjamini-Hochberg FDR control method). The p-values and q-values for up-regulation and down-regulation are reported.
}

\author{
Billy Heung Wing Chang
}

\seealso{
\code{\link{permTestLight}}, \code{\link{wilcox.test}}
}

\examples{
# see the vignette for more examples
# this example is adapted from R GSA package (Efron 2007)

### NOT RUN ###
set.seed(100)
x <- matrix(rnorm(1000*20),ncol=20)
rownames(x) <- paste("g",1:1000,sep="")
dd <- sample(1:1000,size=100)

u <- matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20] <- x[dd,11:20]+u
y <- factor(c(rep('Control',10),rep('Experiment',10)))

results <- wilcoxTest(x, y, tests = "unpaired")
head(results)
}
