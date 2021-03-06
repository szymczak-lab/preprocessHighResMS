preprocessHighResMS
================

This R package implements preprocessing methods for mass spectrometry
data generated with direct infusion high resolution technology
(e.g. FT-ICR).

## Installation

<!-- install.packages(c("ranger", "Boruta", "Umpire", "geoR", "MASS")) -->

First install the dependencies from CRAN:

``` r
## install CRAN dependencies
cran.dep = c("enviPat",
             "FTICRMS",
             "plyr")

cran.dep.to.install = setdiff(cran.dep,
                     installed.packages()[, "Package"])

if (length(cran.dep.to.install) > 0) {
    install.packages(cran.dep.to.install)
}
```

Next install the dependencies from Bioconductor:

``` r
## install Bioconductor dependencies
bioc.dep = c("SummarizedExperiment",
             "xcms")

bioc.dep.to.install = setdiff(bioc.dep,
                              installed.packages()[, "Package"])

if (length(bioc.dep.to.install) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install(bioc.dep.to.install)
}
```

Then install the most recent version of preprocessHighResMS from GitHub:

``` r
library(devtools)
install_github("szymczak-lab/preprocessHighResMS")
```

## Usage

A detailed user guide is available as vignette.
