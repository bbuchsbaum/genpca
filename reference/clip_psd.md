# Clip a symmetric matrix to the PSD cone

Spectral clip: eigen-decompose and set negative eigenvalues to zero.
Unlike
[`ensure_spd()`](https://bbuchsbaum.github.io/genpca/reference/ensure_spd.md)
(a diagonal ridge shift), this preserves the non-negative part of the
spectrum exactly. Requires a dense eigendecomposition, so large sparse
matrices are refused.

## Usage

``` r
clip_psd(M, tol = 1e-06, dense_maxn = 2000L)
```

## Arguments

- M:

  numeric matrix or Matrix::Matrix

- tol:

  tolerance passed to
  [`is_spd()`](https://bbuchsbaum.github.io/genpca/reference/is_spd.md)
  for the fast-path check

- dense_maxn:

  refuse sparse input larger than this (clip densifies)

## Value

a dense Matrix, symmetric PSD
