
<!-- README.md is generated from README.Rmd. Please edit that file -->

# didhetero: Treatment Effect Heterogeneity in Staggered Difference-in-Differences

<!-- badges: start -->
<!-- badges: end -->

The **didhetero** package provides tools to construct doubly robust
uniform confidence bands (UCB) for the group-time conditional average
treatment effect (CATT) function given a pre-treatment covariate of
interest in the staggered difference-in-differences setup of Callaway
and Sant’Anna (2021). The uniform inference procedures are developed by
[Imai, Qin, and Yanagi (2023) “Doubly Robust Uniform Confidence Bands
for Group-Time Conditional Average Treatment Effects in
Difference-in-Differences”](https://doi.org/10.48550/arXiv.2305.02185).

## Installation

Get the package from GitHub:

``` r
# install.packages("devtools") # if needed
devtools::install_github("tkhdyanagi/didhetero", build_vignettes = TRUE)
```

## Vignette

For more details, see the package vignette with:

``` r
# Getting Started with the didhetero Package
vignette("didhetero")
```

## References

- Callaway, B., & Sant’Anna, P. H. (2021). Difference-in-differences
  with multiple time periods. Journal of Econometrics, 225(2), 200-230.
  [Link](https://doi.org/10.1016/j.jeconom.2020.12.001)

- Imai, S., Qin, L., & Yanagi, T. (2023). Doubly Robust Uniform
  Confidence Bands for Group-Time Conditional Average Treatment Effects
  in Difference-in-Differences. arXiv preprint arXiv:2305.02185.
  [Link](https://doi.org/10.48550/arXiv.2305.02185)
