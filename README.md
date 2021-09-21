
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HTEBand

<!-- badges: start -->

<!-- badges: end -->

The goal of HTEBand is to construct uniform confidence bands for various
types of heterogeneous treatment effect parameters.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("remotes") # if not installed
remotes::install_github("koohyun-kwon/HTEBand")
```

## Example

Nonparametric regression via `NpregBand()`:

``` r
library(HTEBand)
library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.4     ✓ dplyr   1.0.7
#> ✓ tidyr   1.1.3     ✓ stringr 1.4.0
#> ✓ readr   2.0.1     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
x <- seq(-1, 1, length.out = 500)
y <- x^2 + rnorm(500, 0, 1/4)
cb.res <- NpregBand(y, x, 2, 0.95, "L", n.eval = 25, c.method = "supp")
#> Residual calculation... Done 
#> Optimal bandwidth calculation... Done 
#> CB construction... 
#> 1 / 25 
#> 2 / 25 
#> 3 / 25 
#> 4 / 25 
#> 5 / 25 
#> 6 / 25 
#> 7 / 25 
#> 8 / 25 
#> 9 / 25 
#> 10 / 25 
#> 11 / 25 
#> 12 / 25 
#> 13 / 25 
#> 14 / 25 
#> 15 / 25 
#> 16 / 25 
#> 17 / 25 
#> 18 / 25 
#> 19 / 25 
#> 20 / 25 
#> 21 / 25 
#> 22 / 25 
#> 23 / 25 
#> 24 / 25 
#> 25 / 25 
#> Done
cb.res$fx <- (cb.res$eval)^2
ggplot(data = cb.res) + geom_line(aes(x = eval, y = cb.lower)) +
  geom_line(aes(x = eval, y = cb.upper)) + geom_line(aes(x = eval, y = fx), color = "red")
```

<img src="man/figures/README-example-1.png" width="100%" />
