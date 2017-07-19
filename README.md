[![Build Status](https://travis-ci.org/billyhw/GSALightning.svg?branch=master)](https://travis-ci.org/billyhw/GSALightning)
[![Build status](https://ci.appveyor.com/api/projects/status/dtokc1ms18i0c020/branch/master?svg=true)](https://ci.appveyor.com/project/billyhw/gsalightning/branch/master)
[![codecov](https://codecov.io/gh/billyhw/GSALightning/branch/master/graph/badge.svg)](https://codecov.io/gh/billyhw/GSALightning)

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
devtools::install_github("billyhw/GSALightning")
```

Since the package contains some large example gene expression data sets, the download time may take 1-3 minutes. Please kindly be patient.

If the R vignette is also desired, use instead:

```{r}
library(devtools)
install_github("billyhw/GSALightning", build_vignette = TRUE)
```

Note the vignette will take some time to build (4-5 minutes) due to some large-scale analysis examples included in the vignette.

The Bioconductor stable release of GSA-Lightning is at:

https://bioconductor.org/packages/release/bioc/html/GSALightning.html

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
