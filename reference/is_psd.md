# Test positive semi-definiteness (relative tolerance)

`is_psd()` is TRUE when `A` is symmetric and every eigenvalue exceeds
`-rtol * max(abs(diag(A)))`; `is_pd()` is TRUE when every eigenvalue
exceeds `+rtol * max(abs(diag(A)))`. Both are shifted Cholesky probes,
so large sparse matrices never need an eigendecomposition. `is_spd()` is
a deprecated alias of `is_psd()` kept for internal callers (its `tol` is
the relative tolerance).

## Usage

``` r
is_psd(A, rtol = .metric_rtol_default())

is_pd(A, rtol = .metric_rtol_default())

is_spd(A, tol = .metric_rtol_default())
```

## Arguments

- A:

  numeric matrix or Matrix::Matrix

- rtol:

  relative tolerance (default `sqrt(.Machine$double.eps)`)

- tol:

  relative tolerance (deprecated name; same as `rtol`)

## Value

logical
