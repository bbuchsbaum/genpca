# Fast generalized matrix decomposition (dense/sparse dispatch)

Computes the generalized SVD of X with row metric Q and column metric R,
equivalent to the eigendecomposition used by
[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) with
`method = "spectra"`. Uses primal (p \<= n) or dual (n \< p) formulation
to minimize computation, with optional Cholesky caching for repeated
calls.

## Usage

``` r
gmd_fast_cpp(
  X,
  Q,
  R,
  k,
  tol = 1e-09,
  maxit = 1000L,
  seed = 1234L,
  topk = TRUE,
  cache = TRUE
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

  tolerance for filtering near-zero singular values. Default 1e-9.

- maxit:

  maximum iterations (ignored, kept for API compatibility)

- seed:

  random seed (ignored, kept for API compatibility)

- topk:

  logical; use top-k symmetric eigen via ARPACK when available. Defaults
  to TRUE. Set to FALSE to force full eigendecomposition.

- cache:

  logical; cache Cholesky factors across calls when constraints are
  dense. Defaults to TRUE. Use
  [`gmd_clear_cache`](https://bbuchsbaum.github.io/genpca/reference/gmd_clear_cache.md)
  to clear.

## Value

A list with components:

- u:

  n x k matrix of scores (left singular vectors scaled)

- v:

  p x k matrix of components (right singular vectors, R-scaled)

- d:

  length-k vector of singular values

- k:

  number of components returned (may be \< requested if rank-deficient)

## When is this fast

This implementation is faster than the "eigen" method when:

- `k << min(n, p)`: Only top-k eigenvalues needed (uses ARPACK)

- Repeated calls with same Q or R: Cholesky factors are cached

- Large matrices where full eigendecomposition is expensive

For small matrices or when `k` is close to `min(n, p)`, the overhead of
iterative methods may make "eigen" faster.

## See also

[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) for
the high-level interface,
[`gmd_clear_cache`](https://bbuchsbaum.github.io/genpca/reference/gmd_clear_cache.md)
to clear Cholesky cache
