The ‘ratioOfQsprays’ package
================
Stéphane Laurent
2024-04-11

*Fractions of multivariate polynomials with rational coefficients.*

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/ratioOfQsprays/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/ratioOfQsprays/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

------------------------------------------------------------------------

The **qspray** package allows arithmetic (and more) on multivariate
polynomials with rational coefficients. Based on this one, the
**ratioOfQsprays** package allows to manipulate fractions of
multivariate polynomials with rational coefficients.

## Creating a `ratioOfQsprays`

A `ratioOfQsprays` object represents a fraction of two multivariate
polynomial with rational coefficients. Such polynomials are represented
by `qspray` objects. The easiest way to create a `ratioOfQsprays` is to
introduce the variables of the polynomials with the `qlone` function,
and then to build a `qspray` numerator and a `qspray` denominator with
arithmetic operations. For example:

``` r
library(ratioOfQsprays)
f <- function(x1, x2, x3) {
  (2*x1^2 + x2*x3) / (4*x1 - 3*x3 + 1)
}
# variables:
x1 <- qlone(1)
x2 <- qlone(2)
x3 <- qlone(3)
# the 'ratioOfQsprays':
( roq <- f(x1, x2, x3) )
## [ 2*x^2 + yz ]  %//%  [ 4*x - 3*z + 1 ]
```

Arithmetic on `ratioOfQsprays` objects is available:

``` r
roq^2
## [ 4*x^4 + 4*x^2yz + y^2z^2 ]  %//%  [ 16*x^2 - 24*xz + 8*x + 9*z^2 - 6*z + 1 ]
roq - roq
## [ 0 ]  %//%  [ 16*x^2 - 24*xz + 8*x + 9*z^2 - 6*z + 1 ]
1 / roq
## [ 4*x - 3*z + 1 ]  %//%  [ 2*x^2 + yz ]
roq + (x2 + x3)/x1
## [ 2*x^3 + xyz + 4*xy + 4*xz - 3*yz + y - 3*z^2 + z ]  %//%  [ 4*x^2 - 3*xz + x ]
```

Rational numbers and polynomials are coercable to `ratioOfQsprays`
objects, and you can also perform arithmetic operations between a
`ratioOfQsprays` and such an object:

``` r
2 * roq
## [ 4*x^2 + 2*yz ]  %//%  [ 4*x - 3*z + 1 ]
"1/2" * roq
## [ x^2 + 1/2*yz ]  %//%  [ 4*x - 3*z + 1 ]
roq + gmp::as.bigq("7/3") 
## [ 2*x^2 + 28/3*x + yz - 7*z + 7/3 ]  %//%  [ 4*x - 3*z + 1 ]
x1 + roq + x3
## [ 6*x^2 + xz + x + yz - 3*z^2 + z ]  %//%  [ 4*x - 3*z + 1 ]
```

## Evaluating a `ratioOfQsprays`

Use `evalRatioOfQsprays` to evaluate a `ratioOfQsprays`. This function
returns a `bigq` number:

``` r
library(gmp) # rational numbers
x <- c("4", "3", "2/5")
evalRatioOfQsprays(roq, x)
## Big Rational ('bigq') :
## [1] 166/79
x <- as.bigq(x)
evalRatioOfQsprays(roq, x)
## Big Rational ('bigq') :
## [1] 166/79
f(x[1], x[2], x[3])
## Big Rational ('bigq') :
## [1] 166/79
```
