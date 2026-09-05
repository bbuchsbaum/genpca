# Get (and cache) a *lower* Cholesky factor for a dense SPD matrix

Get (and cache) a *lower* Cholesky factor for a dense SPD matrix

## Usage

``` r
get_chol_lower_dense(A)
```

## Arguments

- A:

  numeric or dense Matrix (SPD). Sparse input is an error (it is
  factored sparsely elsewhere), never densified here.

## Value

a base numeric matrix L (lower triangular) with A = L %\*% t(L)
