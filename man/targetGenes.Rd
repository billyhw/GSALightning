\name{targetGenes}
\alias{targetGenes}
\docType{data}
\title{
Target Genes of Distal Regulatory Elements
}
\description{
A list containing the target genes of 104636 distal regulatory elements from the human genome. The original list is available in the supplementary data of the reference stated below. The gene names has been tranformed to gene symbol already.
}
\usage{data("targetGenes")}
\value{
A list of distal regulatory elements and their target genes.
}
\source{
http://nar.oxfordjournals.org/content/early/2013/09/03/nar.gkt785/suppl/DC1
}
\references{

Lu, Yulan, Yuanpeng Zhou, and Weidong Tian. "Combining Hi-C data with phylogenetic correlation to predict the target genes of distal regulatory elements in human genome." Nucleic acids research (2013): gkt785.

}
\examples{
data(targetGenes)
str(targetGenes)
}
