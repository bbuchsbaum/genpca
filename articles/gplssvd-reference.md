# Generalized PLS-SVD: Explicit Whitening Reference

This vignette provides a minimal, explicit-whitening reference
implementation for the Generalized PLS-SVD (GPLSSVD) formulation,
matching the operator used in
[`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md).
It is intended for validating small problems and for contributors to
cross-check the mathematical identities.

## Quick practical use

If you just need to run generalized PLS on two blocks, you do not have
to touch the reference code below. The high-level call is:

``` r
set.seed(123)
X <- matrix(rnorm(150*8), 150, 8)
Y <- matrix(rnorm(150*5), 150, 5)
row_wt <- diag(runif(150, 0.5, 1.5))
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

Use this when you want a working analysis; use the reference below when
you want to verify identities.

## Notation

The GPLSSVD decomposes the relationship between two data blocks X (N x
I) and Y (N x J) with optional row and column metrics:

- **MX, MY**: Row metrics (N x N) – weight observations differently for
  X and Y blocks
- **WX, WY**: Column metrics (I x I and J x J) – encode variable
  relationships within each block
- **p, q**: Generalized singular vectors (saliences) satisfying p’WX p =
  I, q’WY q = I
- **Fi, Fj**: Factor scores (component loadings scaled by singular
  values)
- **Lx, Ly**: Latent variables (projections of data onto components)
- **d**: Singular values of the whitened cross-product matrix

``` r
suppressPackageStartupMessages(library(Matrix))

make_sqrt_mats <- function(W, n) {
  if (is.null(W)) {
    list(sqrt = Diagonal(n), invsqrt = Diagonal(n), full = Diagonal(n))
  } else if (is.numeric(W) && length(W) == n) {
    d  <- pmax(as.numeric(W), 0)
    ds <- sqrt(d)
    list(sqrt = Diagonal(x = ds),
         invsqrt = Diagonal(x = ifelse(ds > 0, 1/ds, 0)),
         full = Diagonal(x = d))
  } else if (inherits(W, "diagonalMatrix")) {
    d  <- pmax(as.numeric(diag(W)), 0)
    ds <- sqrt(d)
    list(sqrt = Diagonal(x = ds),
         invsqrt = Diagonal(x = ifelse(ds > 0, 1/ds, 0)),
         full = W)
  } else {
    Wd <- as.matrix(forceSymmetric(Matrix(W)))
    es <- eigen(Wd, symmetric = TRUE)
    Q  <- es$vectors
    lam<- pmax(es$values, 0)
    S  <- Q %*% diag(sqrt(lam), nrow = length(lam)) %*% t(Q)
    IS <- Q %*% diag(ifelse(lam > 0, 1/sqrt(lam), 0), nrow = length(lam)) %*% t(Q)
    list(sqrt = Matrix(S, sparse = FALSE),
         invsqrt = Matrix(IS, sparse = FALSE),
         full = Matrix(Wd, sparse = FALSE))
  }
}

dense_gplssvd_ref <- function(X, Y,
                              MX = NULL, MY = NULL,
                              WX = NULL, WY = NULL,
                              k = NULL,
                              center = FALSE, scale = FALSE) {
  N <- nrow(X); I <- ncol(X); J <- ncol(Y)
  stopifnot(nrow(Y) == N)

  cs <- function(A, do_center, do_scale) {
    A <- as.matrix(A)
    if (is.null(dim(A))) A <- matrix(A, nrow = length(A), ncol = 1)
    cen <- if (isTRUE(do_center)) colMeans(A) else rep(0, ncol(A))
    if (isTRUE(do_center)) A <- A - matrix(rep(cen, each = nrow(A)), nrow(A))
    if (isTRUE(do_scale)) {
      s <- sqrt(colSums(A^2) / pmax(nrow(A) - 1, 1))
      s[s == 0] <- 1
      A <- A %*% diag(1/as.numeric(s), ncol(A))
    } else s <- rep(1, ncol(A))
    list(A = A, center = cen, scale = s)
  }

  Xcs <- cs(X, center, scale); X <- Xcs$A
  Ycs <- cs(Y, center, scale); Y <- Ycs$A

  MXm <- make_sqrt_mats(MX, N)
  MYm <- make_sqrt_mats(MY, N)
  WXm <- make_sqrt_mats(WX, I)
  WYm <- make_sqrt_mats(WY, J)

  Xe <- MXm$sqrt %*% X %*% WXm$sqrt
  Ye <- MYm$sqrt %*% Y %*% WYm$sqrt
  S  <- crossprod(Xe, Ye)

  sv <- svd(as.matrix(S))
  if (!is.null(k)) {
    k <- min(k, length(sv$d))
    sv$u <- sv$u[, seq_len(k), drop = FALSE]
    sv$v <- sv$v[, seq_len(k), drop = FALSE]
    sv$d <- sv$d[seq_len(k)]
  }

  p <- WXm$invsqrt %*% sv$u
  q <- WYm$invsqrt %*% sv$v
  Fi <- WXm$full %*% (p %*% diag(sv$d, nrow = length(sv$d)))
  Fj <- WYm$full %*% (q %*% diag(sv$d, nrow = length(sv$d)))
  Lx <- MXm$sqrt %*% (X %*% (WXm$full %*% p))
  Ly <- MYm$sqrt %*% (Y %*% (WYm$full %*% q))

  list(d = sv$d, u = sv$u, v = sv$v, p = p, q = q, fi = Fi, fj = Fj, lx = Lx, ly = Ly,
       S = S,
       center = list(X = Xcs$center, Y = Ycs$center),
       scale  = list(X = Xcs$scale,  Y = Ycs$scale))
}
```

Example usage:

``` r
set.seed(1)
N <- 20; I <- 8; J <- 6
X <- matrix(rnorm(N*I), N, I)
Y <- matrix(rnorm(N*J), N, J)
MX <- diag(runif(N, .5, 1.5))
MY <- diag(runif(N, .5, 1.5))
WX <- diag(runif(I, .5, 1.5))
WY <- diag(runif(J, .5, 1.5))

ref <- dense_gplssvd_ref(X, Y, MX, MY, WX, WY, k = 3, center = TRUE, scale = FALSE)

op  <- genpca::gplssvd_op(X, Y, XLW = MX, YLW = MY, XRW = WX, YRW = WY,
                          k = 3, center = TRUE, scale = FALSE)

cat("Singular values match:\n")
#> Singular values match:
print(all.equal(ref$d, op$d, tolerance = 1e-6))
#> [1] TRUE

cat("\nLatent variable cross-products equal singular values:\n")
#> 
#> Latent variable cross-products equal singular values:
print(all.equal(diag(crossprod(op$lx, op$ly)), op$d, tolerance = 1e-6))
#> [1] TRUE

cat("\nSingular values:\n")
#> 
#> Singular values:
print(round(op$d, 4))
#> [1] 22.0777 19.9684 12.8428
```

This reference is intentionally small and dense to highlight the exact
whitening relationships; the package’s operator-based implementation
avoids materializing large matrices for scalability and sparsity.

## Where next

See
[`vignette("gpca-basics", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-basics.md)
for a getting-started guide and
[`vignette("gpca-metrics", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-metrics.md)
for metric recipes.
