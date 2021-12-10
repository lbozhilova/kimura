
<!-- README.md is generated from README.Rmd. Please edit that file -->

\*\* NOTE: This is an untested and poorly documented toy
reimplementation of
[lbozhilova/Kimura-Distribution)](https://github.com/lbozhilova/Kimura-Distribution).
\*\*

# kimura

<!-- badges: start -->
<!-- badges: end -->

The goal of `kimura` is to fit Kimura dsitributions to heteroplasmy data
and test for evidence of selection pressure.

## Installation

You can install the development version of `kimura` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
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
#> D = 0.12798, p = 0.43400, b = 0.74170, p-value = 0.992
#> alternative hypothesis: one-sided
```

There is also some additional functionality for generating data from the
Kimura distirbution.

``` r
# Initialise Kimura parameters
p <- 0.3250
b <- 0.01257

# Probability of allele loss
dkimura(0, p, b)
#> [1] 0.6667266
# Probability of fixing an allele
dkimura(1, p, b)
#> [1] 0.3167281

# Kimura(p, d) CDF at 0.1 intervals
pkimura(seq(0, 1, 0.1), p, b)
#>  [1] 0.6667266 0.6683807 0.6700356 0.6716903 0.6733450 0.6749996 0.6766540
#>  [8] 0.6783084 0.6799627 0.6816169 1.0000000

# Random number generation
rkimura(10, p, b)
#>  [1] 0 1 1 1 0 0 1 0 0 1
```

## References

Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels. “The
distribution of mitochondrial DNA heteroplasmy due to random genetic
drift.” The American Journal of Human Genetics 83.5 (2008): 582-593.

Kimura, Motoo. “Solution of a process of random genetic drift with a
continuous model.” Proceedings of the National Academy of Sciences of
the United States of America 41.3 (1955): 144.
