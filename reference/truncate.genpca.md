# Truncate a genpca fit to fewer components

Returns a new `genpca` object retaining only the first `ncomp`
components. All component-indexed slots (`v`, `s`, `sdev`, `ov`, `ou`,
`u`, `propv`, `cumv`) are sliced consistently; the preprocessing object
and constraint matrices are carried over unchanged.

## Usage

``` r
# S3 method for class 'genpca'
truncate(x, ncomp)
```

## Arguments

- x:

  A `genpca` object.

- ncomp:

  Number of components to retain (a positive integer no larger than
  `ncomp(x)`).

## Value

A `genpca` object with `ncomp` components.

## See also

[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md),
[`reconstruct.genpca()`](https://bbuchsbaum.github.io/genpca/reference/reconstruct.genpca.md)

## Examples

``` r
X <- matrix(rnorm(60), 15, 4)
fit <- genpca(X, ncomp = 4)
fit2 <- truncate(fit, 2)
multivarious::ncomp(fit2)
#> [1] 2
```
