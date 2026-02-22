# Get (and cache) a \*lower\* Cholesky factor for a dense SPD matrix

Get (and cache) a \*lower\* Cholesky factor for a dense SPD matrix

## Usage

``` r
get_chol_lower_dense(A)
```

## Arguments

- A:

  numeric or dense Matrix (SPD). If sparse, falls back to dense.

## Value

a base numeric matrix L (lower triangular) with A = L
