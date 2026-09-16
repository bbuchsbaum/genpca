# Reconstruct data from a genpca fit

Reconstructs (an approximation of) the original data from a `genpca` fit
as `ou[, comp] %*% diag(d[comp]) %*% t(ov[, comp])`, followed by the
inverse of the preprocessing transform. With all components and full
rank this recovers the original data.

## Usage

``` r
# S3 method for class 'genpca'
reconstruct(
  x,
  comp = 1:multivarious::ncomp(x),
  rowind = NULL,
  colind = NULL,
  ...
)
```

## Arguments

- x:

  A `genpca` object.

- comp:

  Integer vector of components to use (default: all).

- rowind:

  Optional integer vector of rows to reconstruct (default: all).

- colind:

  Optional integer vector of columns to reconstruct (default: all). The
  inverse preprocessing transform is applied to the selected columns.

- ...:

  Ignored.

## Value

A numeric matrix of dimension `length(rowind) x length(colind)`.

## See also

[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md),
[`truncate.genpca()`](https://bbuchsbaum.github.io/genpca/reference/truncate.genpca.md)

## Examples

``` r
X <- matrix(rnorm(60), 15, 4)
fit <- genpca(X, ncomp = 4, preproc = multivarious::center())
max(abs(reconstruct(fit) - X)) # ~ 0 at full rank
#> [1] 8.881784e-16
```
