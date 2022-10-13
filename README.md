
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kimura

<!-- badges: start -->
<!-- badges: end -->

The goal of `kimura` is to fit Kimura distributions to heteroplasmy data
and test for evidence of selection pressure.

## Installation

You can install the development version of `kimura` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library("devtools")
devtools::install_github("lbozhilova/kimura")
```

## Example

The principal purpose of this package is to perform hypothesis test for
the presence of selection pressure in mtDNA heteroplasmy distributions
(Wonnapinij et al., 2008). This can be done like so:

``` r
library(kimura)

# Load some heteroplasmy data
h <- c(0.06, 0.08, 0.27, 0.37, 0.40, 0.45, 0.56, 0.61, 0.75, 0.79)

# Carry out test for selection 
test_kimura(h)
#> 
#>  Monte Carlo Kolmogorov-Smirnov
#> 
#> data:  h and Kimura(0.434, 0.7417)
#> D = 0.12798, p = 0.43400, b = 0.74170, p-value = 0.989
#> alternative hypothesis: one-sided
```

There is also some additional functionality for generating data from the
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

## References

Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels. “The
distribution of mitochondrial DNA heteroplasmy due to random genetic
drift.” The American Journal of Human Genetics 83.5 (2008): 582-593.

Kimura, Motoo. “Solution of a process of random genetic drift with a
continuous model.” Proceedings of the National Academy of Sciences of
the United States of America 41.3 (1955): 144.

This package started out as a toy reimplementation of
[lbozhilova/Kimura-Distribution](https://github.com/lbozhilova/Kimura-Distribution).
