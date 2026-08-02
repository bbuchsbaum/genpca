# Generalized PLS via Implicit Operator (PLS-SVD / GPLSSVD)

Canonical (two-block) generalized PLS using sparse-friendly implicit
matrix-vector products. Solves the SVD of the operator \\S = Xe' Ye\\
without materializing \\Xe = Mx^{1/2} X Ax^{1/2}\\ or \\Ye = My^{1/2} Y
Ay^{1/2}\\.

## Usage

``` r
genpls(
  X,
  Y,
  Ax = NULL,
  Ay = NULL,
  Mx = NULL,
  My = NULL,
  ncomp = 2,
  preproc_x = multivarious::pass(),
  preproc_y = multivarious::pass(),
  svd_backend = c("RSpectra", "irlba"),
  svd_opts = list(tol = 1e-07, maxitr = 1000),
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric or Matrix, n x p.

- Y:

  Numeric or Matrix, n x q. Must have same n as `X`.

- Ax:

  Column metric for X (W_X): vector/diagonal/matrix; `NULL` means
  identity.

- Ay:

  Column metric for Y (W_Y): vector/diagonal/matrix; `NULL` means
  identity.

- Mx:

  Row metric for X (M_X): vector/diagonal/matrix; `NULL` means identity.

- My:

  Row metric for Y (M_Y): vector/diagonal/matrix; `NULL` means identity.

- ncomp:

  Number of components to extract (rank-k). Default 2.

- preproc_x, preproc_y:

  Optional `multivarious` preprocessors (e.g., `center()`). Defaults to
  [`multivarious::pass()`](https://bbuchsbaum.github.io/multivarious/reference/pass.html)
  (no-op).

- svd_backend:

  Character, one of `"RSpectra"` (default) or `"irlba"` for the
  iterative SVD. This choice only matters for larger problems: whenever
  both `X` and `Y` have at most 64 columns after preprocessing, the
  operator materializes `S` densely and computes a direct
  [`svd()`](https://rdrr.io/r/base/svd.html), ignoring `svd_backend`
  entirely (see
  [`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md)).

- svd_opts:

  List of options passed to the SVD backend, e.g., `tol`, `maxitr`.

- verbose:

  Logical; print brief progress messages.

## Value

An object of class `c("genpls", "cross_projector", "projector")` with:

- vx, vy:

  X- and Y- projection weights (stored in cross_projector) such that
  `project(fit, X) = X \%*\% vx` and
  `project(fit, Y, source = "Y") = Y \%*\% vy` recover the latent
  variables in the ambient (non-whitened) metric. Algebraically
  `vx = W_X p = fi \%*\% diag(1/d)` and
  `vy = W_Y q = fj \%*\% diag(1/d)` (see Details).

- d:

  singular values of \\S = Xe' Ye\\ (attached field)

- p, q:

  generalized weights \\W_X^{-1/2} u\\, \\W_Y^{-1/2} v\\ (attached)

- fi, fj:

  variable/component scores \\W_X p D\\, \\W_Y q D\\ (attached)

- lx, ly:

  row latent variables \\M_X^{1/2} X W_X p\\, \\M_Y^{1/2} Y W_Y q\\
  (attached)

- metrics:

  the supplied metrics (attached)

- ncomp:

  Number of components actually extracted. The underlying operator may
  return fewer than the requested `ncomp` (e.g. when `ncomp` exceeds
  `min(ncol(X), ncol(Y))`); this field reflects the actual count, not
  the request.

- backend:

  The `svd_backend` value passed in (for reference only; see the
  `svd_backend` argument for when it is actually used).

- preproc_x, preproc_y:

  The fitted `multivarious` preprocessing objects for `X` and `Y`,
  stored on the `cross_projector`.

## Details

This follows the GPLSSVD/PLS-SVD formulation (Beaton, eqs. 10–14): the
top `ncomp` singular triplets of S are computed by iterative SVD on the
linear maps v -\> S v and u -\> S^T u, implemented with metric Cholesky
multiplies/solves when possible. Works with dense or sparse `Matrix`
inputs and constraint metrics.

Returns a
[`multivarious::cross_projector`](https://bbuchsbaum.github.io/multivarious/reference/cross_projector.html)
with X-/Y-weights (vx, vy) chosen to provide natural projection of new
data (`X %*% vx`, `Y %*% vy`). Additional GPLSSVD quantities are
attached to the object for access: singular values `d`, generalized
weights `p`, `q`, variable scores `fi`, `fj`, and row latent variables
`lx`, `ly`.

`genpls()` maximizes the covariance between latent variables of `X` and
`Y` under the (Mx, Ax, My, Ay) metrics by computing the SVD of \\S = Xe'
Ye\\, where \\Xe = Mx^{1/2} X Ax^{1/2}\\ and \\Ye = My^{1/2} Y
Ay^{1/2}\\.

`project()` on a fitted object returns latent variables in the ambient
(original data) metric, i.e. `X \%*\% vx = X W_X p`. This differs from
the attached `lx = Mx^{1/2} X W_X p`, which lives in the row-whitened
metric, by the factor \\Mx^{1/2}\\: `lx` and `project(fit, X)` are equal
only when `Mx = I`. New-row projection necessarily uses `project()`'s
ambient-metric convention, because a training-row metric `Mx` has no
natural extension to out-of-sample rows.

**Metric naming.** `genpls()`'s row/column metric arguments (`Mx`, `My`
for rows; `Ax`, `Ay` for columns) follow the same M/A convention as
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md).
Internally they are forwarded to
[`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md),
which uses the names `XLW`/`YLW` (left/row weights, i.e. `Mx`/`My`) and
`XRW`/`YRW` (right/column weights, i.e. `Ax`/`Ay`).

## References

Beaton, D. (2020). Generalized eigen, singular value, and partial least
squares decompositions: The GSVD package. (Eqs. 10-14).
arXiv:2010.14734.

## Examples

``` r
if (requireNamespace("RSpectra", quietly = TRUE) &&
    requireNamespace("multivarious", quietly = TRUE)) {
  set.seed(1)
  n <- 100; p <- 40; q <- 30
  X <- matrix(rnorm(n*p), n, p)
  Y <- matrix(rnorm(n*q), n, q)
  w <- runif(n); w <- w/sum(w)
  Mx <- My <- Matrix::Diagonal(x = w)
  fit <- genpls(X, Y, Mx = Mx, My = My, ncomp = 2,
                preproc_x = multivarious::center(),
                preproc_y = multivarious::center())
  fit$d  # singular values
}
#> [1] 1.401768 1.340951
```
