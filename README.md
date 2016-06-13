## Synopsis

GSA-Lightning provides a fast implementation of permutation-based gene set
analysis for two-sample problem. This package is particularly useful when testing
simultaneously a large number of gene sets, or when a large number of permutations
is necessary for more accurate p-values estimation.

## Motivation

The computational speed of many gene set analysis methods can be slow due to the
computationally demanding permutation step. GSA-Lightning is a fast implementation
of permutation-based gene set analysis. GSA-Lightning achieves significant speedup compared to existing
R implementations, particularly when the number of gene sets and permutations are large.

## Installation

We recommend using the R "devtools" package to install GSALightning from Github. To install devtools, in an R prompt type:

```{r}
install.packages("devtools")
```

Then use the following to install GSA-Lightning:

```{r}
library(devtools) 
install_github("billyhw/GSALightning")
```

GSA-Lightning is also available at Bioconductor:

https://bioconductor.org/packages/release/bioc/html/GSALightning.html

## Code

Please consult the R help page of GSALightning for examples:

```{r}
? GSALight
```

## Contributors

Billy Heung Wing Chang

## Reference

Billy Heung Wing Chang and Weidong Tian (2016) GSA-Lightning: Ultra-fast Permutation-based Gene Set Analysis. Bioinformatics. doi: 10.1093/bioinformatics/btw349

http://bioinformatics.oxfordjournals.org/content/early/2016/06/11/bioinformatics.btw349.abstract

## License

GPL (>=2)
