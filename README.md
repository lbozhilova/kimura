
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kimura

<!-- badges: start -->
<!-- badges: end -->

The goal of `kimura` is to provide functionality related to the Kimura distribution from statistical genetics (Kimura, 1955), often applied to the study of mitochondrial DNA heteroplasmy.

## Installation

You can install the development version of `kimura` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library("devtools")
devtools::install_github("lbozhilova/kimura")
```

## Example

This package contains functionality for extracting probabilities and generating data from the
Kimura distribution.

``` r
# Initialise Kimura parameters
p <- 0.6
b <- 0.95

# Probability of allele loss
dkimura(0, p, b)
#> [1] 1.989732e-05
# Probability of fixing an allele
dkimura(1, p, b)
#> [1] 9.694172e-06

# Kimura(p, d) CDF at 0.1 intervals
pkimura(seq(0, 1, 0.1), p, b)
#>  [1] 1.989732e-05 2.165073e-05 1.345668e-04 3.560559e-03 3.773170e-02
#>  [6] 1.840206e-01 4.908528e-01 8.136989e-01 9.711726e-01 9.989777e-01
#> [11] 1.000000e+00

# Random number generation
rkimura(10, p, b)
#>  [1] 0.6837089 0.7312714 0.5824325 0.4693752 0.3754452 0.5464930 0.6789625
#>  [8] 0.6312104 0.7449587 0.5605141
```


This package also contains functionality for fitting a Kimura distribution
to a heteroplasmy sample using the method of moments. This has in the past
been used to carry out 
a Kolmogorov-Smirnov test to assess discrepancy between the sample and 
the fitted distribution (Wonnapinij et al., 2008). 

❗ However, since (Wonnapinij et al., 2008), it has been found that a 
method-of-moments fit does not provide the best-fitting version of the Kimura
distribution for a given heteroplasmy sample. Such a fit will often 
produce a distribution that is a bad fit to the data, and therefore a 
*false positive signal of selection*, whereas a differently fitted 
distribution would given no such signal. This issue is described in
(Giannakis et al., 2023).

❗ It is *not recommended* to assess selection by fitting a Kimura distribution to
snapshot heteroplasmy data -- time series data is more appropriate for this 
task. This repository https://github.com/StochasticBiology/heteroplasmy-analysis/
contains functions to identify better Kimura fits than method of moments.


## References

Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels. “The
distribution of mitochondrial DNA heteroplasmy due to random genetic
drift.” The American Journal of Human Genetics 83.5 (2008): 582-593.

Kimura, Motoo. “Solution of a process of random genetic drift with a
continuous model.” Proceedings of the National Academy of Sciences of
the United States of America 41.3 (1955): 144.

Giannakis, K., Broz, A.K., Sloan, D.B. and Johnston, I.G., 2023. Avoiding misleading estimates using mtDNA heteroplasmy statistics to study bottleneck size and selection. G3: Genes, Genomes, Genetics, 13(6), p.jkad068.

This package started out as a toy reimplementation of
[lbozhilova/Kimura-Distribution](https://github.com/lbozhilova/Kimura-Distribution).
