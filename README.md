The ‘ratioOfQsprays’ package
================
Stéphane Laurent
2024-04-13

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
introduce the variables of the polynomials with the `qlone` function
(from the **qspray** package), and then to build a `qspray` numerator
and a `qspray` denominator with the arithmetic operations. For example:

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
## [ 2*x^2 + y.z ]  %//%  [ 4*x - 3*z + 1 ]
```

Arithmetic on `ratioOfQsprays` objects is available:

``` r
roq^2
## [ 4*x^4 + 4*x^2.y.z + y^2.z^2 ]  %//%  [ 16*x^2 - 24*x.z + 8*x + 9*z^2 - 6*z + 1 ]
roq - roq
## [ 0 ]
1 / roq
## [ 4*x - 3*z + 1 ]  %//%  [ 2*x^2 + y.z ]
2*roq + (x2 + x3)/x1
## [ 4*x^3 + 2*x.y.z + 4*x.y + 4*x.z - 3*y.z + y - 3*z^2 + z ]  %//%  [ 4*x^2 - 3*x.z + x ]
```

Rational numbers and `qspray` polynomials are coercable to
`ratioOfQsprays` objects, and then you can also perform arithmetic
operations between a `ratioOfQsprays` and such an object:

``` r
2 * roq
## [ 4*x^2 + 2*y.z ]  %//%  [ 4*x - 3*z + 1 ]
"1/2" * roq
## [ x^2 + 1/2*y.z ]  %//%  [ 4*x - 3*z + 1 ]
roq + gmp::as.bigq("7/3") 
## [ 2*x^2 + 28/3*x + y.z - 7*z + 7/3 ]  %//%  [ 4*x - 3*z + 1 ]
x1 + roq + x3^2
## [ 6*x^2 + 4*x.z^2 - 3*x.z + x + y.z - 3*z^3 + z^2 ]  %//%  [ 4*x - 3*z + 1 ]
```

The result of an arithmetic operation is always an irreducible fraction.
To perform this step, the C++ library **CGAL** is used to compute a
greatest common divisor of the numerator and the denominator of the
possibly non-reduced fraction resulting from the arithmetic operation,
and then to divide both of them by this greatest common divisor. This is
very efficient in general.

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

It is also possible to substitute for a subset of the variables, with
the help of the function `substituteRatioOfQsprays`. You have to
indicate the variables you don’t want to replace with `NA`:

``` r
x <- c(NA, "3", "2/5")
substituteRatioOfQsprays(roq, x)
## [ 2*x^2 + 6/5 ]  %//%  [ 4*x - 1/5 ]
x <- as.bigq(x)
f(x1, x[2], x[3])
## [ 2*x^2 + 6/5 ]  %//%  [ 4*x - 1/5 ]
```

And it is possible to convert a `ratioOfQsprays` to a function which is
evaluated by **Ryacas**:

``` r
fyac <- as.function(roq)
fyac("4", "3", "2/5") # = evalRatioOfQsprays(roq, c("4", "3", "2/5"))
## [1] "166/79"
```

Actually you can pass some literal variables to this function:

``` r
fyac("x", "3", "2/5") # = substituteRatioOfQsprays(roq, c(NA, "3", "2/5"))
## [1] "(2*(5*x^2+3))/(20*x-1)"
fyac("x", "y", "z")   # = roq
## [1] "(y*z+2*x^2)/(4*x-3*z+1)"
fyac("x", "x", "x")
## [1] "(3*x^2)/(x+1)"
```

Complex numbers and allowed; the imaginary unit is denoted by `I`. See
the **Yacas** documentation for more information.

``` r
fyac("Sqrt(2)", "2 + 2*I", "3")
## [1] "Complex(10/(Sqrt(32)-8),6/(Sqrt(32)-8))"
```

You can get numerical approximations by setting the option `N=TRUE` in
`as.function`:

``` r
fyacN <- as.function(roq, N = TRUE)
fyacN("4", "3", "2/5") 
## [1] 2.101266
fyacN("x", "3", "2/5")
## expression((2 * (5 * x^2 + 3))/(20 * x - 1))
fyacN("Sqrt(2)", "2 + 2*I", "3")
## [1] -4.267767-2.56066i
```

## Querying a `ratioOfQsprays`

A couple of functions to query a `ratioOfQsprays` are available:

``` r
getNumerator(roq)
## 2*x^2 + y.z
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

Now, if you perform an arithmetic operation between `roq` *at first
position* and an another `ratioOfQsprays` or an object coercable to a
`ratioOfQsprays`, these show options are passed to the result if
possible:

``` r
roq + (x1 + 1)/x2
## [ 2*a1^2.a2 + 4*a1^2 - 3*a1.a3 + 5*a1 + a2^2.a3 - 3*a3 + 1 ] / [ 4*a1.a2 - 3*a2.a3 + a2 ]
```

Here it is always possible to pass the options to the result. An obvious
example for which this is not always possible is when you use three
letters for the variables, e.g.

``` r
showRatioOfQspraysOption(roq, "showQspray") <- showQsprayXYZ(c("A", "B", "C"))
roq
## [ 2*A^2 + B.C ] / [ 4*A - 3*C + 1 ]
```

but then you add a `ratioOfQsprays` containing the fourth variable:

``` r
roq + x4/(x4 + 1)
## [ 2*A1^2.A4 + 2*A1^2 + 4*A1.A4 + A2.A3.A4 + A2.A3 - 3*A3.A4 + A4 ] / [ 4*A1.A4 + 4*A1 - 3*A3.A4 - 3*A3 + A4 + 1 ]
```

Obviously it is not possible to denote the resulting fraction of
polynomials with the letters `A`, `B` and `C`. The solution I adopted
consists in taking the first of these letters and to index it. The same
method is used for the `qspray` polynomials.
