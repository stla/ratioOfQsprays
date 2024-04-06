The ‘ratioOfQsprays’ package
================
Stéphane Laurent
2024-04-07

*Fractions of multivariate polynomials with rational coefficients.*

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/ratioOfQsprays/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/ratioOfQsprays/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

------------------------------------------------------------------------

The **qspray** package allows arithmetic (and more) on multivariate
polynomials with rational coefficients. Based on this one, the
**ratioOfQsprays** package allows to manipulate fractions of
multivariate polynomials with rational coefficients.

A `symbolicQspray` object represents a multivariate polynomial whose
coefficients are fractions of polynomials. For example:

``` r
library(ratioOfQsprays)
f <- function(x1, x2, x3) {
  2*x1^2 / (3*x3 + 1)  +  (x2 * x3) / (4*x1 + x2)
}
# variables:
x1 <- qlone(1)
x2 <- qlone(2)
x3 <- qlone(3)
# the 'ratioOfQsprays':
( roq <- f(x1, x2, x3) )
## [ 8*x1^3 + 2*x1^2.x2 + 3*x2.x3^2 + x2.x3 ]  %//%  [ 12*x1.x3 + 4*x1 + 3*x2.x3 + x2 ]
```

Arithmetic on `ratioOfQsprays` objects is available:

``` r
roq^2
## [ 64*x1^6 + 32*x1^5.x2 + 4*x1^4.x2^2 + 48*x1^3.x2.x3^2 + 16*x1^3.x2.x3 + 12*x1^2.x2^2.x3^2 + 4*x1^2.x2^2.x3 + 9*x2^2.x3^4 + 6*x2^2.x3^3 + x2^2.x3^2 ]  %//%  [ 144*x1^2.x3^2 + 96*x1^2.x3 + 16*x1^2 + 72*x1.x2.x3^2 + 48*x1.x2.x3 + 8*x1.x2 + 9*x2^2.x3^2 + 6*x2^2.x3 + x2^2 ]
roq - roq
## [0]
1 / roq
## [ 12*x1.x3 + 4*x1 + 3*x2.x3 + x2 ]  %//%  [ 8*x1^3 + 2*x1^2.x2 + 3*x2.x3^2 + x2.x3 ]
```

And evaluation:

``` r
library(gmp) # rational numbers
x <- c(as.bigq(4), as.bigq(3), as.bigq("2/5"))
evalRatioOfQsprays(roq, x)
## Big Rational ('bigq') :
## [1] 15266/1045
f(x[1], x[2], x[3])
## Big Rational ('bigq') :
## [1] 15266/1045
```
