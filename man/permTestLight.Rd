\name{permTestLight}
\alias{permTestLight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Single-Gene Permutation Test
}
\description{
A fast permutation-testing procedure for two-sample single-gene differential expression analysis.
}
\usage{
permTestLight(eset, fac, nperm, tests = c('unpaired','paired'), method = c("mean", "absmean"), 
	      npermBreaks = 2000, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
The expression matrix. Each row is a gene, and each column is a subject/sample. The gene names must be presented as the row names.
}
  \item{fac}{
Subject labels, for unpaired T-test, either a factor or something that can be coerced into a factor (e.g. 0 and 1, Experiment and Control). For paired T-test,  fac must be an integer vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions.
}
  \item{nperm}{
Number of permutations. If unspecified, nperm will be set as the total number of gene sets divided by 0.05 times 2. This should be sufficient to estimate accurate p-values for Bonferroni Correction under significance level alpha = 0.05.
}
  \item{tests}{
The tests to performed. Can be either the default "unpaired" for unpaired T-tests or "paired" for paired T-tests.
}
  \item{method}{
Modification of the T-test statistics for each individual gene for hypothesis testing. The default "mean" option uses the typical T-statistics. The other option "absmean" uses the absolute value of the T-statistics (which results in a two-sided test).
}
  \item{npermBreaks}{
The batch size. When the number of permutation nperm is large, the permutations are broken into batches, and that permutation are performed with the batches sequentially run. Default is 2000.
}
  \item{verbose}{
Should the progress be reported? Default = TRUE.
}
}
\details{
The speed performance is sensitive to npermBreaks. Setting npermBreaks small can save memory, but will take longer to run. Setting npermBreaks large can speed up the process, but may run into memory issues. If the function is running slow, consider increase npermBreaks. If the function is running into memory issues, consider reducing npermBreaks. The default 2000 typically can provide a reasonable balance between speed and memory.
}
\value{
A data frame with the p-values, the q-values (via Benjamini-Hochberg FDR control method), the gene statistics, and the gene set size. If method = 'mean', then the p-values and q-values for up-regulation and down-regulation are reported. Method = 'absmean' corresponds to a two-sided test of absolute changes in expression, hence only one set of p-values and one set of q-values will be reported.
}
\references{
BHW Chang and W Tian (2015). GSA-Lightning: Ultra Fast Permutation-based Gene Set Analysis. (Submitted)
}
\author{
Billy Heung Wing Chang
}

\seealso{
\code{\link{GSALight}}
}
\examples{

# see the vignette for more examples
# this example is adapted from R GSA package (Efron 2007)

set.seed(100)
x <- matrix(rnorm(1000*20),ncol=20)
rownames(x) <- paste("g",1:1000,sep="")
dd <- sample(1:1000,size=100)

u <- matrix(2*rnorm(100),ncol=10,nrow=100)
x[dd,11:20] <- x[dd,11:20]+u
y <- factor(c(rep('Control',10),rep('Experiment',10)))

results <- permTestLight(x, y, nperm = 1000, method = 'mean')
head(results)

}
