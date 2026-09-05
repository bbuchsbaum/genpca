# Clip a symmetric matrix to the PSD cone

Spectral clip: eigen-decompose and set negative eigenvalues to zero.
Unlike
[`ensure_spd()`](https://bbuchsbaum.github.io/genpca/reference/ensure_spd.md)
(a diagonal ridge shift), this preserves the non-negative part of the
spectrum exactly. The output has no negative eigenvalue beyond
reconstruction roundoff: the only fast path is an exact (unshifted)
Cholesky success, which proves positive definiteness. Requires a dense
eigendecomposition, so large sparse matrices are refused.

## Usage

``` r
clip_psd(M, tol = NULL, dense_maxn = 2000L, name = "M")
```

## Arguments

- M:

  numeric matrix or Matrix::Matrix

- tol:

  unused; kept for call compatibility

- dense_maxn:

  refuse sparse input larger than this (clip densifies)

- name:

  label used in error messages

## Value

a symmetric Matrix, PSD
