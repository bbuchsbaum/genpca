# GPCA at Scale and Special Cases

## Backend selection

- `method = "eigen"`: small/medium dense problems; forms full matrices.
- `method = "spectra"`: large or sparse; matrix-free iterations
  (requires C++ build).
- `method = "deflation"`: need only a few components; minimal memory.

## Sparse workflow (spectra)

``` r
set.seed(42)
X_sparse <- rsparsematrix(800, 300, density = 0.005)
fit_sp <- genpca(X_sparse, ncomp = 5, method = "spectra",
                 preproc = multivarious::pass(),
                 constraints_remedy = "ridge")
fit_sp$sdev
```

## Covariance-only GPCA

When you already have the cross-product matrix C = X’MX,
[`genpca_cov()`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
avoids working with the full data:

``` r
set.seed(123)
n <- 100; p <- 15
X <- matrix(rnorm(n * p), n, p)
M <- diag(runif(n, 0.8, 1.2))
A <- diag(runif(p, 0.7, 1.3))
C <- t(X) %*% M %*% X
fit_cov <- genpca_cov(C, R = A, ncomp = 5, method = "gmd")
fit_cov$d
#> [1] 13.80217 12.42550 11.92054 11.15895 10.96560
```

## Out-of-sample projection

Fit on training rows, then project held-out observations into the same
space:

``` r
set.seed(7)
X <- matrix(rnorm(200 * 30), 200, 30)
fit <- genpca(X[1:150, ], ncomp = 4, preproc = multivarious::center())
scores_test <- multivarious::project(fit, X[151:200, ])
head(scores_test)
#>             PC1        PC2        PC3        PC4
#> [1,]  1.9323426 -0.4080526 -0.1924407 -0.6673048
#> [2,]  0.3498745  0.6490485  0.2579339  0.9996621
#> [3,]  0.9590113  0.9305126 -1.4603670 -1.1496702
#> [4,]  0.1125973  1.0845757 -0.2493420 -1.2509707
#> [5,] -1.5977372 -0.9139256  1.5647260  0.1034587
#> [6,] -0.2072342 -0.7558773  1.2817970 -1.0316951
```

## Performance tips

Keep data sparse when possible; avoid centering dense copies if you plan
to use `spectra`. Use `constraints_remedy = "ridge"` for empirical
metrics. For big problems, limit `ncomp` to what you truly need and
consider the covariance route if `n` is huge but `p` is moderate.

## Where next

See
[`vignette("gpca-metrics", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-metrics.md)
for building metrics and
[`vignette("gpca-basics", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-basics.md)
for quick-start examples.
