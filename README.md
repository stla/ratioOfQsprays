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
## [ 0 ]
1 / roq
## [ 4*x - 3*z + 1 ]  %//%  [ 2*x^2 + yz ]
2*roq + (x2 + x3)/x1
## [ 4*x^3 + 2*xyz + 4*xy + 4*xz - 3*yz + y - 3*z^2 + z ]  %//%  [ 4*x^2 - 3*xz + x ]
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

## Querying a `ratioOfQsprays`

A couple of functions to query a `ratioOfQsprays` are available:

``` r
getNumerator(roq)
## 2*x^2 + yz
getDenominator(roq)
## 4*x - 3*z + 1
numberOfVariables(roq)
## [1] 3
isConstant(roq)
## [1] FALSE
isConstant(roq / roq)
## [1] TRUE
isUnivariate(roq)
## [1] FALSE
isUnivariate(x1 / (x1^2 + 1))
## [1] TRUE
isPolynomial(roq)
## [1] FALSE
isPolynomial((x1^2 - x2^2)/(x1 - x2))
## [1] TRUE
```

## Showing a `ratioOfQsprays`

As you have seen, the variables of `roq` are denoted by `x`, `y`, `z`.
This is the default of showing a `ratioOfQsprays` which have no more
than three variables. If it has more than three variables, the variables
are denoted by `x1`, `x2`, `x3`, …:

``` r
x4 <- qlone(4)
roq / x4
## [ 2*x1^2 + x2.x3 ]  %//%  [ 4*x1.x4 - 3*x3.x4 + x4 ]
```

It is possible to control the way a `ratioOfQsprays` is printed. For
example, let’s say you want to print `roq` by using `a1`, `a2`, `a3` for
the variables and you want to change the symbol for the quotient:

``` r
showRatioOfQspraysOption(roq, "x") <- "a"
showRatioOfQspraysOption(roq, "quotientBar") <- " / " 
roq
## [ 2*a1^2 + a2.a3 ] / [ 4*a1 - 3*a3 + 1 ]
```

Now, if you perform an arithmetic operation between `roq` at first
position and an another `ratioOfQsprays` or an object coercable to a
`ratioOfQsprays`, these show options are preserved:

``` r
roq + (x1 + 1)/x2
## [ 2*a1^2.a2 + 4*a1^2 - 3*a1.a3 + 5*a1 + a2^2.a3 - 3*a3 + 1 ] / [ 4*a1.a2 - 3*a2.a3 + a2 ]
```
