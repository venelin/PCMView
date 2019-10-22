
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCMView

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/venelin/PCMView.svg?branch=master)](https://travis-ci.com/venelin/PCMView)
<!-- badges: end -->

The goal of PCMView is to provide visualization functions for various
objects generated with the R-packages PCMBase and PCMFit.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("venelin/PCMView")
```

## Example

Here is how to plot the search path from a recurisive clade partition
run:

``` r
library(PCMView)
#> Registered S3 method overwritten by 'treeio':
#>   method     from
#>   root.phylo ape

lstPlots <- PlotRCPHistory(
   fitBMWithShiftsToMammalData, layout = "fan")
cowplot::plot_grid(plotlist = lstPlots[!sapply(lstPlots, is.null)])
```

<img src="man/figures/README-PlotRCPHistory-1.png" width="100%" />
