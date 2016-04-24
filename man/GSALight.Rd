\name{GSALight}
\alias{GSALight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Permutation-based Gene Set Analysis
}
\description{
GSALight is a fast implementation of two-sample permutation-based gene set analysis. It supports the mean or absolute mean version of the permutation T-test as implemeneted in the GSA package. Restandardization is also supported.
}
\usage{
GSALight(eset, fac, gs, nperm = NULL, tests = c('unpaired','paired'), method = c("mean", "absmean"), 
	 minsize = 1, maxsize = Inf, restandardize = TRUE,
	 npermBreaks = 2000, rmGSGenes = c('stop', 'gene', 'gs'), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
The expression matrix. Each row is a gene, and each column is a subject/sample. The gene names must be presented as the row names.
}
  \item{fac}{
Subject labels, for unpaired T-test, either a factor or something that can be coerced into a factor (e.g. 0 and 1, Experiment and Control). For paired T-test,  fac must be an integer vector of 1,-1,2,-2,..., where each number represents a pair, and the sign represents the conditions.
}
  \item{gs}{
The gene sets. Most commonly, gs is a list, where each element (named after the gene set) is a vector of genes belonging to the gene set. It can also be a binary sparse matrix, where each row is a gene set, and each column is a gene. For each row (i.e. a gene set), the row entry is 1 if the corresponding gene belongs to the gene set, and 0 otherwise. Alternatively it can be a data table, where the first column (must be named "geneSet") contains the names of the gene set, and the second column (must be named "gene") are the gene set genes. 
}
  \item{nperm}{
Number of permutations. If unspecified, nperm will be set as the total number of gene sets divided by 0.05 times 2. This should be sufficient to estimate accurate p-values for Bonferroni Correction under significance level alpha = 0.05.
}
  \item{tests}{
The tests to performed. Can be either the default "unpaired" for unpaired T-tests or "paired" for paired T-tests.
}
  \item{method}{
The method to combine the T-test statistics for each individual gene into a gene set test statistics. The default "maxmean" uses the maxmean statistics proposed by Efron (2007). Other options are "mean", i.e. the mean of the statistics of the genes inside a gene set, and "absmean", the mean of the absolute value of the statistics.
}
  \item{minsize}{
Minimum gene set size (default = 1, i.e. even gene set with a single gene is allowed).
}
  \item{maxsize}{
Maximum gene set size (default = Inf, i.e. no upper limit).
}
  \item{restandardize}{
Should restandardization be performed? This is typically recommended to avoid excessive number of significant gene sets.
}
  \item{npermBreaks}{
The batch size. When the number of permutation nperm is large, the permutations are broken into batches, and that permutation are performed with the batches sequentially run. Default is 2000.
}
  \item{rmGSGenes}{
What to do if there gene sets with genes without expression measurements. Default is 'stop', and will return an error if there are gene set genes with missing expression measurements. If rmGSGenes is "gene", genes with missing expression measurements are removed from the gene sets. If rmGSGenes = 'gs', gene sets with non-measured genes are removed.
}
  \item{verbose}{
Should the progress be reported? Default = TRUE.
}
}
\details{
The speed performance of GSALight is sensitive to npermBreaks. Setting npermBreaks small can save memory, but will take longer to run. Setting npermBreaks large can speed up GSALight, but may run into memory issues. If GSALight is running slow, consider increasing npermBreaks. If GSALight is running into memory issues, consider reducing npermBreaks. The default 2000 can typically provide a reasonable balance between speed and memory.
}
\value{
A data frame with the p-values, the q-values (via Benjamini-Hochberg FDR control method), the gene set statistics, and the gene set size. If method = 'mean', then the p-values and q-values for up-regulation and down-regulation are reported. Method = 'absmean' corresponds to a one-sided test of absolute changes in expression, hence only one set of p-values and one set of q-values will be reported.
}
\references{
BHW Chang and W Tian (2015). GSA-Lightning: Ultra Fast Permutation-based Gene Set Analysis. (Submitted)

B Efron and RJ Tibshirani (2007). "On testing the significance of sets of genes." The annals of applied statistics 1(1):107-129
}
\author{
Billy Heung Wing Chang
}

\seealso{
\code{\link{permTestLight}}, \code{\link{targetGenes}}
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

#create some random gene sets
genesets=vector("list",50)
for(i in 1:50){
 genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
}
names(genesets)=paste("set",as.character(1:50),sep="")

GSAmean <- GSALight(x, y, genesets, nperm = 1000, method = 'mean', restandardize = FALSE, rmGSGenes = 'gene')
GSAabs <- GSALight(x, y, genesets, nperm = 1000, method = 'absmean', restandardize = FALSE, rmGSGenes = 'gene')

head(GSAmean)
head(GSAabs)

}
