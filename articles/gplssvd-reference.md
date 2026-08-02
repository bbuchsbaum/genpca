# Generalized PLS-SVD: Explicit Whitening Reference

This vignette has two audiences. End users who just want to run
generalized PLS on two blocks should read the **Quick practical use**
section and stop. Contributors who want to verify the whitening
identities behind
[`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md)
should continue to the explicit reference implementation below.

## Quick practical use

``` r

set.seed(123)
N <- 150
X <- matrix(rnorm(N * 8), N, 8)
Y <- matrix(rnorm(N * 5), N, 5)
row_wt   <- diag(runif(N, 0.5, 1.5))
col_wt_x <- diag(runif(8, 0.8, 1.2))
col_wt_y <- diag(runif(5, 0.8, 1.2))

fit <- genpls(X, Y, ncomp = 2,
              preproc_x = multivarious::center(),
              preproc_y = multivarious::center(),
              Mx = row_wt, My = row_wt,
              Ax = col_wt_x, Ay = col_wt_y)
round(fit$d, 3)
#> [1] 57.572 47.290
```

![Singular values of the generalized
cross-product.](gplssvd-reference_files/figure-html/quick-plot-1.png)

Singular values of the generalized cross-product.

That is enough to fit and inspect a model. The remainder of this
vignette is for contributors verifying the math.

## Notation

GPLSSVD decomposes the relationship between two data blocks `X`
(`N x I`) and `Y` (`N x J`) with optional row and column metrics:

- `MX`, `MY`: row metrics (`N x N`) – weight observations differently
  for the two blocks
- `WX`, `WY`: column metrics (`I x I` and `J x J`) – encode within-block
  variable relationships
- `p`, `q`: generalized singular vectors (saliences) satisfying
  `p' WX p = I`, `q' WY q = I`
- `Fi`, `Fj`: factor scores (loadings scaled by singular values)
- `Lx`, `Ly`: latent variables (data projections onto components)
- `d`: singular values of the whitened cross-product matrix

## Reference implementation

The function below builds the whitened cross-product
`S = (M_X^{1/2} X W_X^{1/2})' (M_Y^{1/2} Y W_Y^{1/2})` explicitly and
runs a dense SVD. It is deliberately the most literal transcription of
the algebra above, with no attention to speed or to input types beyond
what the check below needs; the package’s operator-based path avoids
materialising the whitened matrices at all.

Everything rests on one helper. Each metric enters through its symmetric
PSD square root, and singular metrics need the *pseudo*-inverse of that
root – zero eigenvalues stay zero rather than blowing up:

``` r

psd_sqrt <- function(W, n) {
  if (is.null(W)) return(list(h = diag(n), hinv = diag(n), full = diag(n)))
  W    <- as.matrix(W)
  e    <- eigen(W, symmetric = TRUE)
  lam  <- pmax(e$values, 0)                    # clip round-off negatives
  half <- function(f) e$vectors %*% (f * t(e$vectors))
  list(h    = half(sqrt(lam)),                 # W^{1/2}
       hinv = half(ifelse(lam > 0, 1 / sqrt(lam), 0)),  # W^{-1/2}, pseudo
       full = W)
}
```

The decomposition itself is then a whitening, one SVD, and an
unwhitening:

``` r

dense_gplssvd_ref <- function(X, Y, MX = NULL, MY = NULL,
                              WX = NULL, WY = NULL, k = NULL,
                              center = FALSE, scale = FALSE) {
  X <- scale(as.matrix(X), center = center, scale = scale)
  Y <- scale(as.matrix(Y), center = center, scale = scale)
  stopifnot(nrow(X) == nrow(Y))

  mx <- psd_sqrt(MX, nrow(X)); wx <- psd_sqrt(WX, ncol(X))
  my <- psd_sqrt(MY, nrow(Y)); wy <- psd_sqrt(WY, ncol(Y))

  # whiten both blocks, then SVD their cross-product
  S  <- crossprod(mx$h %*% X %*% wx$h, my$h %*% Y %*% wy$h)
  sv <- svd(S)
  keep <- seq_len(if (is.null(k)) length(sv$d) else min(k, length(sv$d)))

  # unwhiten: saliences are W^{-1/2} u, so that p' WX p = I
  p <- wx$hinv %*% sv$u[, keep, drop = FALSE]
  q <- wy$hinv %*% sv$v[, keep, drop = FALSE]
  D <- diag(sv$d[keep], nrow = length(keep))

  list(d  = sv$d[keep], p = p, q = q,
       fi = wx$full %*% p %*% D,             # factor scores
       fj = wy$full %*% q %*% D,
       lx = mx$h %*% X %*% wx$full %*% p,    # latent variables
       ly = my$h %*% Y %*% wy$full %*% q)
}
```

## Cross-checking the operator

Run the reference on a small block, run
[`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md)
with the same metrics, and compare:

``` r

set.seed(1)
N <- 20; I <- 8; J <- 6
X  <- matrix(rnorm(N * I), N, I)
Y  <- matrix(rnorm(N * J), N, J)
MX <- diag(runif(N, .5, 1.5))
MY <- diag(runif(N, .5, 1.5))
WX <- diag(runif(I, .5, 1.5))
WY <- diag(runif(J, .5, 1.5))

ref <- dense_gplssvd_ref(X, Y, MX, MY, WX, WY,
                         k = 3, center = TRUE, scale = FALSE)
op  <- gplssvd_op(X, Y,
                  XLW = MX, YLW = MY,
                  XRW = WX, YRW = WY,
                  k = 3, center = TRUE, scale = FALSE)

all.equal(ref$d, op$d, tolerance = 1e-6)
#> [1] TRUE
all.equal(diag(crossprod(op$lx, op$ly)), op$d, tolerance = 1e-6)
#> [1] TRUE
round(op$d, 4)
#> [1] 22.0777 19.9684 12.8428
```

![Reference vs operator singular values agree to plotting precision
(left). Latent variables show the expected diagonal cross-product
structure
(right).](gplssvd-reference_files/figure-html/ref-vs-op-plot-1.png)

Reference vs operator singular values agree to plotting precision
(left). Latent variables show the expected diagonal cross-product
structure (right).

The diagonal of `t(Lx) %*% Ly` recovers the singular values, as the
GPLSSVD identity guarantees.

## Where next

See
[`vignette("genpca")`](https://bbuchsbaum.github.io/genpca/articles/genpca.md)
for a getting-started walkthrough and
[`vignette("gpca-metrics")`](https://bbuchsbaum.github.io/genpca/articles/gpca-metrics.md)
for metric recipes that apply to both GPCA and GPLSSVD.
