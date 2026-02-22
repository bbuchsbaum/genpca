# Generalized PLS-SVD via Implicit Operator (memory-safe)

Compute the top-k singular triplets of \\S = Xe' Ye\\ without
materializing the whitened matrices \\Xe = Mx^{1/2} X Wx^{1/2}\\, \\Ye =
My^{1/2} Y Wy^{1/2}\\. Works with dense/sparse constraints.

## Usage

``` r
gplssvd_op(
  X,
  Y,
  XLW = NULL,
  YLW = NULL,
  XRW = NULL,
  YRW = NULL,
  k = 2,
  center = FALSE,
  scale = FALSE,
  svd_backend = c("RSpectra", "irlba"),
  svd_opts = list(tol = 1e-07, maxitr = 1000)
)
```

## Arguments

- X:

  n x I matrix (numeric or Matrix)

- Y:

  n x J matrix (numeric or Matrix)

- XLW:

  Row metric for X (M_X): NULL/identity, numeric length-n,
  diagonalMatrix, or PSD Matrix

- YLW:

  Row metric for Y (M_Y)

- XRW:

  Column metric for X (W_X)

- YRW:

  Column metric for Y (W_Y)

- k:

  Number of components

- center, scale:

  Logical; pre-center/scale columns of X, Y before metrics

- svd_backend:

  One of "RSpectra" (default) or "irlba"

- svd_opts:

  List of options for the backend (e.g., tol, maxitr)

## Value

list with elements d, u, v, p, q, fi, fj, lx, ly, k, dims, center, scale
