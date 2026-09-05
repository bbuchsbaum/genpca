# Generalized matrix decomposition via partial SVD of the whitened operator

Computes the generalized SVD of X with row metric Q and column metric R,
equivalent to the eigendecomposition used by
[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) with
`method = "eigen"`. The metrics are factored once (`Q = F_Q F_Q'`,
`R = F_R F_R'`; diagonal, dense or sparse Cholesky, or an eigen factor
for singular metrics) and the top-k singular triplets of the implicit
operator `F_Q' X F_R` are computed with eigencore; a dense SVD is used
when few components are not requested or the iterative solver does not
converge. `gmd_fast_cpp()` is an alias kept for existing callers.

## Usage

``` r
gmd_spectra(
  X,
  Q,
  R,
  k,
  tol = 1e-09,
  maxit = 1000L,
  seed = 1234L,
  topk = TRUE,
  cache = TRUE,
  auto_topk = TRUE,
  topk_ratio = 0.08,
  topk_min_dim = 200L,
  diag_fast = TRUE,
  rank_rtol = 1e-06,
  metric_rtol = .metric_rtol_default(),
  dense_maxn = 5000L
)

gmd_fast_cpp(
  X,
  Q,
  R,
  k,
  tol = 1e-09,
  maxit = 1000L,
  seed = 1234L,
  topk = TRUE,
  cache = TRUE,
  auto_topk = TRUE,
  topk_ratio = 0.08,
  topk_min_dim = 200L,
  diag_fast = TRUE,
  rank_rtol = 1e-06,
  metric_rtol = .metric_rtol_default(),
  dense_maxn = 5000L
)
```

## Arguments

- X:

  numeric matrix (n x p)

- Q, R:

  constraints (weights/metrics) for rows/cols. Must be symmetric
  positive (semi-)definite. Can be dense matrices, sparse matrices, or
  diagonal matrices.

- k:

  number of components to extract (must be \>= 1 and \<= min(n, p))

- tol:

  convergence tolerance of the iterative solver. Default 1e-9.

- maxit:

  unused (kept for compatibility).

- seed:

  unused (kept for compatibility); results do not depend on the R random
  stream.

- topk:

  logical; use the iterative top-k solver when `k < min(n, p)`. Set to
  FALSE to force a dense SVD of the whitened operator.

- cache:

  logical; cache dense Cholesky factors across calls. Defaults to TRUE.
  Use
  [`gmd_clear_cache`](https://bbuchsbaum.github.io/genpca/reference/gmd_clear_cache.md)
  to clear.

- auto_topk:

  logical; when TRUE (default), use top-k only when `k/min(n,p)` is
  small and `min(n,p)` is large enough.

- topk_ratio:

  threshold used by `auto_topk`. If `k/min(n,p) <= topk_ratio`, top-k is
  used. Default 0.08.

- topk_min_dim:

  minimum `min(n,p)` required before top-k is used under `auto_topk`.
  Default 200.

- diag_fast:

  logical; if TRUE (default) and both constraints are diagonal, use a
  weighted-SVD fast path.

- rank_rtol:

  relative cutoff on singular values: components with
  `d_j <= rank_rtol * d_1` are dropped. Default 1e-6.

- metric_rtol:

  relative tolerance for metric validation and null-space detection.
  Default `sqrt(.Machine$double.eps)`.

- dense_maxn:

  a singular general metric on the small side of X needs a dense
  eigendecomposition; refuse it above this many rows (the `maxeig`
  argument of
  [`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)).
  A singular metric on the large side is never factored: the solver
  switches to the symmetric small-side formulation, in which that metric
  only appears in products.

## Value

A list with components:

- u:

  n x k matrix of metric-weighted scores `Q ou D`

- v:

  p x k matrix of components `R ov`

- ou,ov:

  metric-orthonormal factors

- d:

  length-k vector of singular values

- k:

  number of components returned (may be \< requested if rank-deficient)

## When is this fast

- `k << min(n, p)`: only the top-k triplets are computed

- Repeated calls with the same dense Q or R: Cholesky factors are cached

- Sparse metrics: only sparse factors and products are formed

A positive definite metric on the big side of X costs one Cholesky of
that dimension; a singular one is never factored (symmetric small-side
form).

## See also

[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) for
the high-level interface,
[`gmd_clear_cache`](https://bbuchsbaum.github.io/genpca/reference/gmd_clear_cache.md)
to clear the Cholesky cache
