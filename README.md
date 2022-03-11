
<!-- README.md is generated from README.Rmd. Please edit that file -->

\*\* NOTE: This is an untested and poorly documented toy
reimplementation of
[lbozhilova/Kimura-Distribution](https://github.com/lbozhilova/Kimura-Distribution).\*\*

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
#> D = 0.12798, p = 0.43400, b = 0.74170, p-value = 0.991
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
#>  [1] 0.7022604 0.6391005 0.5597187 0.5749873 0.7201949 0.6307230 0.5723443
#>  [8] 0.6263145 0.7428140 0.7362410
```

## References

Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels. “The
distribution of mitochondrial DNA heteroplasmy due to random genetic
drift.” The American Journal of Human Genetics 83.5 (2008): 582-593.

Kimura, Motoo. “Solution of a process of random genetic drift with a
continuous model.” Proceedings of the National Academy of Sciences of
the United States of America 41.3 (1955): 144.
