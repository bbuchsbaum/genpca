# Ensure SPD (sparse-friendly)

Force a symmetric matrix to be symmetric positive definite (SPD). Uses a
Gershgorin-based diagonal shift; falls back to nearPD for small dense
matrices.

## Usage

``` r
ensure_spd(M, tol = 1e-06, nearpd_maxn = 2000L)
```

## Arguments

- M:

  numeric matrix or Matrix::Matrix

- tol:

  jitter tolerance (default 1e-8)

- nearpd_maxn:

  only use nearPD when n \<= nearpd_maxn and matrix is dense

## Value

a Matrix object (sparse stays sparse when possible)
