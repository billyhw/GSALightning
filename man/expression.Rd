\name{expression}
\alias{expression}
\docType{data}
\title{
Breast Cancer Data from The Cancer Genome Atlas (TCGA)
}
\description{

This is a breast cancer expression data set obtained from The Cancer Genome Atlas (TCGA) consortium and processed by the Pan-Cancer project group. The data was downloaded originally using the ELMER Bioconductor package. The gene names have been converted to gene symbols in this data.

}
\usage{data("expression")}
\format{
A matrix with 20299 genes and 1218 subjects, of which 114 are controls.
}

\references{
The Cancer Genome Atlas (2012) Comprehensive molecular characterization of human colon and rectal cancer. Nature 487:330-337.

Weinstein, John N., et al. "The cancer genome atlas pan-cancer analysis project." Nature genetics 45.10 (2013): 1113-1120.

Yao, L., et al. "Inferring regulatory element landscapes and transcription factor networks from cancer methylomes." Genome biology 16.1 (2015): 105-105.
}
\examples{
data(expression)
str(expression)
}

