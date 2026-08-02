# Generalized PLS-SVD via Implicit Operator (memory-safe)

Compute the top-k singular triplets of \\S = Xe' Ye\\ without
materializing the whitened matrices \\Xe = Mx^{1/2} X Wx^{1/2}\\, \\Ye =
My^{1/2} Y Wy^{1/2}\\ when doing so would densify sparse data. When the
whitening is sparsity-preserving (identity/diagonal metrics) or the data
are dense, the whitened blocks are precomputed once so each
matrix-vector product in the iterative SVD costs two multiplies.

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

  Number of components. If `k` exceeds `min(ncol(X), ncol(Y))`, a
  warning is issued and `k` is silently truncated to that maximum.

- center, scale:

  Logical; pre-center/scale columns of X, Y before metrics

- svd_backend:

  One of "RSpectra" (default) or "irlba". Ignored whenever both
  `ncol(X) <= 64` and `ncol(Y) <= 64`, in which case `S` is materialized
  densely and solved with
  [`base::svd()`](https://rdrr.io/r/base/svd.html).

- svd_opts:

  List of options for the backend (e.g., tol, maxitr)

## Value

A list with elements:

- d:

  Length-`k` numeric vector of singular values of \\S = Xe' Ye\\.

- u:

  `I x k` matrix; left singular vectors of `S` (orthonormal in the
  Euclidean metric).

- v:

  `J x k` matrix; right singular vectors of `S` (orthonormal in the
  Euclidean metric).

- p:

  `I x k` matrix of generalized X-weights, \\p = W_X^{-1/2} u\\.

- q:

  `J x k` matrix of generalized Y-weights, \\q = W_Y^{-1/2} v\\.

- fi:

  `I x k` matrix of X-variable scores, \\F_i = W_X p D\\ (columns of `p`
  rescaled by the singular values).

- fj:

  `J x k` matrix of Y-variable scores, \\F_j = W_Y q D\\.

- lx:

  `N x k` matrix of X row latent variables, \\L_x = M_X^{1/2} X W_X p\\.

- ly:

  `N x k` matrix of Y row latent variables, \\L_y = M_Y^{1/2} Y W_Y q\\.

- k:

  Integer; number of components actually returned (may be less than the
  requested `k` if it exceeded `min(I, J)`).

- dims:

  A list `list(N, I, J)` with the row count `N` and column counts
  `I = ncol(X)`, `J = ncol(Y)`.

- center:

  A list `list(X, Y)` of the length-`I` / length-`J` column means
  subtracted from `X`/`Y` (all zero when `center = FALSE`).

- scale:

  A list `list(X, Y)` of the length-`I` / length-`J` column scale
  factors divided out of `X`/`Y` (all one when `scale = FALSE`).

## Details

Naming map: this function names its metrics `XLW`/`YLW` (left/row
weights) and `XRW`/`YRW` (right/column weights); these correspond to
`Mx`/`My` (row metrics) and `Ax`/`Ay` (column metrics) in
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)'s
and
[`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md)'s
M/A convention.

## References

Beaton, D. (2020). Generalized eigen, singular value, and partial least
squares decompositions: The GSVD package. arXiv:2010.14734.

Abdi, H. (2007). Partial least square regression PLS-Regression. In N.
Salkind (Ed.), *Encyclopedia of Measurement and Statistics*. Thousand
Oaks, CA: Sage.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(40 * 6), 40, 6)
Y <- matrix(rnorm(40 * 4), 40, 4)
op <- gplssvd_op(X, Y, k = 2, center = TRUE)
round(op$d, 3)
#> [1] 20.428 14.384
```
