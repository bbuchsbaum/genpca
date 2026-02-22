# Generalized PCA on a covariance matrix

Performs Generalized PCA directly on a pre-computed covariance matrix C
with a single variable-side constraint/metric R. This is useful when you
already have C = X'MX or when X is too large to store but C is
manageable. Supports two methods: "gmd" (Allen et al.'s GMD approach,
default) which exactly matches the two-sided
[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md), and
"geigen" (generalized eigenvalue approach) which solves C v = lambda R
v.

## Usage

``` r
genpca_cov(
  C,
  R = NULL,
  ncomp = NULL,
  method = c("gmd", "geigen"),
  constraints_remedy = c("error", "ridge", "clip", "identity"),
  tol = 1e-08,
  verbose = FALSE
)
```

## Arguments

- C:

  A p x p symmetric positive semi-definite covariance matrix. Typically
  C = X'MX where X is the data matrix and M is a row metric.

- R:

  Variable-side constraint/metric. Can be:

  - NULL: Identity matrix (standard PCA on C)

  - A numeric vector of length p: Interpreted as diagonal weights (must
    be non-negative)

  - A p x p symmetric PSD matrix: General metric/smoothing/structure
    penalties

- ncomp:

  Number of components to return. Default is all positive eigenvalues.

- method:

  Character string specifying the method. One of:

  - "gmd" (default): Allen et al.'s GMD approach via eigen decomposition
    of \\R^{1/2} C R^{1/2}\\

  - "geigen": Generalized eigenvalue approach solving C v = lambda R v

- constraints_remedy:

  How to handle slightly non-PSD inputs (for geigen method). One of:

  - "error": Stop with an error if constraints are not PSD

  - "ridge": Add a small ridge to the diagonal to make PSD

  - "clip": Clip negative eigenvalues to zero

  - "identity": Replace with identity matrix

- tol:

  Numerical tolerance for PSD checks and filtering small eigenvalues.
  Default 1e-8.

- verbose:

  Logical. If TRUE, print progress messages. Default FALSE.

## Value

A list with components:

- v:

  p x k matrix of loadings (R-orthonormal eigenvectors)

- d:

  Singular values (square root of eigenvalues lambda)

- lambda:

  Eigenvalues (variances under the R-metric)

- k:

  Number of components returned

- propv:

  Proportion of variance explained by each component

- cumv:

  Cumulative proportion of variance explained

- R_rank:

  Rank of the constraint matrix R

- method:

  The method used ("gmd" or "geigen")

## Details

**Method Selection Guide:**

Use `method = "gmd"` when:

- You need exact equivalence with
  [`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)`(X, M, A)`

- You're following Allen et al.'s GMD framework

- You want consistent results with the two-sided decomposition

Use `method = "geigen"` when:

- You specifically need the generalized eigenvalue formulation

- You're working with legacy code that expects this approach

- Computational efficiency is critical and R is well-conditioned

**Method "gmd" (default):**

This method implements Allen et al.'s GMD approach and exactly matches
the two-sided genpca when C = X'MX. It computes the eigendecomposition
of \\R^{1/2} C R^{1/2}\\ and maps back with \\V = R^{-1/2} Z\\, ensuring
V'RV = I. The total variance is tr(CR) as in Allen's GPCA (Corollary 5).

**Method "geigen":**

This method solves the generalized eigenproblem C v = lambda R v
directly. While mathematically valid, it solves a different optimization
than Allen's GMD and will not, in general, match the two-sided genpca
unless R = I or special commutation conditions hold.

For exact equivalence with genpca(X, M, A), use method="gmd" with C =
X'MX and R = A.

## References

Allen, G. I., Grosenick, L., & Taylor, J. (2014). A Generalized
Least-Squares Matrix Decomposition. Journal of the American Statistical
Association, 109(505), 145-159.

## See also

[`genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.md) for
the standard two-sided GPCA on data matrices,
[`genpls`](https://bbuchsbaum.github.io/genpca/reference/genpls.md) for
generalized partial least squares

## Examples

``` r
# Example 1: Standard PCA on covariance (no constraint)
C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
fit0 <- genpca_cov(C, R=NULL, ncomp=3)
print(fit0$d[1:3])       # first 3 singular values
#> [1] 2.0562689 0.4926162 0.2796596
print(fit0$propv[1:3])   # variance explained by first 3 components
#> [1] 0.92461872 0.05306648 0.01710261

# Example 2: Demonstrating equivalence with genpca
set.seed(123)
X <- matrix(rnorm(50 * 10), 50, 10)
M_diag <- runif(50, 0.5, 1.5)  # row weights
A_diag <- runif(10, 0.5, 2)    # column weights

# Two-sided GPCA
fit_gpca <- genpca(X, M = M_diag, A = A_diag, ncomp = 5,
                   preproc = multivarious::pass())

# Equivalent covariance-based GPCA
C <- crossprod(X, diag(M_diag) %*% X)  # C = X'MX
fit_cov <- genpca_cov(C, R = A_diag, ncomp = 5, method = "gmd")

# These should match exactly
all.equal(fit_gpca$sdev, fit_cov$d, tolerance = 1e-10)
#> [1] TRUE

# Example 3: Variable weights via a diagonal metric (using iris covariance)
C_iris <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
w <- c(1, 1, 0.5, 2)  # emphasize Sepal.Width less, Petal.Width more
fitW <- genpca_cov(C_iris, R = w, ncomp=3, method="gmd")
print(fitW$d[1:3])
#> [1] 1.7972881 0.4903437 0.3173048

# Example 4: Compare GMD and generalized eigenvalue approaches
fit_gmd <- genpca_cov(C_iris, R = w, ncomp=2, method="gmd")
fit_geigen <- genpca_cov(C_iris, R = w, ncomp=2, method="geigen")
# These will generally differ unless R = I
print(paste("GMD singular values:", paste(round(fit_gmd$d, 3), collapse=", ")))
#> [1] "GMD singular values: 1.797, 0.49"
print(paste("GEigen singular values:", paste(round(fit_geigen$d, 3), collapse=", ")))
#> [1] "GEigen singular values: 2.658, 0.497"
```
