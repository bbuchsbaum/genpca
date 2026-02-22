# Check if constraint matrix is identity or purely diagonal with all != 1

Check if constraint matrix is identity or purely diagonal with all != 1

## Usage

``` r
is_identity_or_diag(M, eps = 1e-15)
```

## Arguments

- M:

  A matrix (often \`dsCMatrix\`) or \`NULL\`.

- eps:

  Numeric tolerance

## Value

TRUE if `M` is a (Matrix-based) diagonal with all or partial diag
