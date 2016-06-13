## Motivation

GSA-Lightning is an R package that provides a fast implementation of permutation-based gene set
analysis for two-sample problem. This package is particularly useful when testing
simultaneously a large number of gene sets, or when a large number of permutations
is necessary for more accurate p-values estimation.

## Installation

GSA-Lightning requires R 3.3.0 or above (https://www.r-project.org/).

We recommend using the R "devtools" package to install GSALightning from Github. 

To install devtools, in an R prompt type:

```{r}
install.packages("devtools")
```

To use devtools to install GSA-Lightning, type:

```{r}
library(devtools) 
install_github("billyhw/GSALightning")
```

GSA-Lightning is also available at Bioconductor:

https://bioconductor.org/packages/release/bioc/html/GSALightning.html

Also consider the development version:

http://bioconductor.org/packages/devel/bioc/html/GSALightning.html

## Code

Please consult the R help page of GSA-Lightning for examples:

```{r}
? GSALight
```

## Contributor

Billy Heung Wing Chang (billyheungwing@gmail.com)

## Reference

Billy Heung Wing Chang and Weidong Tian (2016) GSA-Lightning: Ultra-fast Permutation-based Gene Set Analysis. Bioinformatics. doi: 10.1093/bioinformatics/btw349

http://bioinformatics.oxfordjournals.org/content/early/2016/06/11/bioinformatics.btw349.abstract

## License

GPL (>=2)
