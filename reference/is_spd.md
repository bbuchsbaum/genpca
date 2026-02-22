# Check if matrix is SPD

Check if a matrix is symmetric positive semi-definite using Cholesky
decomposition

## Usage

``` r
is_spd(A, tol = 1e-06)
```

## Arguments

- A:

  numeric matrix or Matrix::Matrix

- tol:

  tolerance for numerical checks (unused but kept for compatibility)

## Value

logical TRUE if SPD, FALSE otherwise
