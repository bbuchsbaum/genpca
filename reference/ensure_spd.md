# Ensure SPD (sparse-friendly)

Force a symmetric matrix to be symmetric positive definite: the result
satisfies `is_pd(., rtol = tol)`. Already-PD input is returned
unchanged; otherwise a Gershgorin-based diagonal shift is applied, with
a [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
fallback for small dense matrices and an escalating jitter as a last
resort.

## Usage

``` r
ensure_spd(M, tol = 1e-06, nearpd_maxn = 2000L, name = "M")
```

## Arguments

- M:

  numeric matrix or Matrix::Matrix

- tol:

  relative positive-definiteness margin (default 1e-6)

- nearpd_maxn:

  only use nearPD when n \<= nearpd_maxn and matrix is dense

- name:

  label used in error messages

## Value

a Matrix object (sparse stays sparse when possible)
