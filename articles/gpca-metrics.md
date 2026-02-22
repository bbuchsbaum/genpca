# GPCA Metrics: Building M and A

This vignette collects practical recipes for row/column metrics, plus
notes on SPD remedies and the experimental
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md).

## Why metrics

Metrics encode weighting and correlation. Row metric `M` changes how
observations are compared; column metric `A` changes how variables are
compared. Both must be symmetric positive definite (SPD) for GPCA.

## Runnable example: heteroscedastic diagonals

The simplest non-trivial metrics are diagonal matrices that down-weight
noisy rows or columns:

``` r
set.seed(42)
n <- 60; p <- 20
X <- matrix(rnorm(n * p), n, p)

col_noise_sd <- runif(p, 0.5, 2)
A <- Diagonal(x = 1 / (col_noise_sd^2))
row_noise_sd <- runif(n, 0.7, 1.3)
M <- Diagonal(x = 1 / (row_noise_sd^2))

fit <- genpca(X, M = M, A = A, ncomp = 3, preproc = multivarious::center())
fit$sdev
#> [1] 14.69818 11.24551 10.96887
```

## Recipes (swap in your data)

Temporal AR(1) rows (smooth in time):

``` r
rho <- 0.7; n <- 200
idx <- 0:(n-1)
Sigma_r <- outer(idx, idx, function(i,j) rho^abs(i-j))
M <- solve(Sigma_r + 1e-3 * diag(n))
```

Spatial kernel (rows or cols):

``` r
coords <- as.matrix(expand.grid(x = 1:10, y = 1:10))
d2 <- as.matrix(dist(coords))^2
ell <- 2
K <- exp(-d2 / (2 * ell^2))
A <- solve(K + 1e-3 * diag(nrow(K)))
```

Graph Laplacian (columns):

``` r
W <- bandSparse(30, k = c(-1,0,1), diagonals = list(rep(0.2,29), rep(1,30), rep(0.2,29)))
D <- Diagonal(x = rowSums(W))
L <- D - W
A <- L + 1e-2 * Diagonal(nrow(L))
```

## Learning metrics with `gpca_mle`

When you suspect structured noise but lack good priors,
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
can learn SPD metrics by alternating GPCA with matrix-normal maximum
likelihood:

``` r
fit_mle <- gpca_mle(X, ncomp = 3, max_iter = 10,
                    lambda = 1e-3, scale_fix = "trace",
                    method = "eigen", verbose = FALSE)
fit_mle$A
fit_mle$M
```

Keep `lambda` non-zero, start with small `ncomp`, and inspect warnings
(they mean a metric was repaired).

## SPD remedies in practice

- `constraints_remedy = "ridge"` (default) adds a diagonal jitter to
  keep metrics SPD.
- `"clip"` truncates tiny negative eigenvalues; `"identity"` falls back
  to identity if a metric misbehaves.
- Scale metrics before use (e.g., divide by mean diagonal) to avoid
  ill-conditioning.
- After fitting, check `range(eigen(A)$values)` on small problems to
  verify conditioning.

## Where next

See
[`vignette("gpca-scale", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-scale.md)
for backend choices, sparse workflows, and covariance-only GPCA.
