# Check if matrix is SPD

Check if a matrix is symmetric positive semi-definite (within `tol`) by
attempting a Cholesky factorization of `A + tol*scale*I`.

## Usage

``` r
is_spd(A, tol = 1e-06)
```

## Arguments

- A:

  numeric matrix or Matrix::Matrix

- tol:

  relative tolerance: eigenvalues above `-tol * max(abs(diag(A)))` are
  treated as non-negative, so PSD-but-singular metrics are accepted

## Value

logical TRUE if symmetric PSD within tolerance, FALSE otherwise
