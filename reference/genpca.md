# Generalised Principal Components Analysis (GPCA)

Implements the Generalised Least‑Squares Matrix Decomposition of Allen,
Grosenick & Taylor (2014) for data observed in a \*\*row\*\*
inner‑product space M and a \*\*column\*\* inner‑product space A.
Setting M = I_n, A = I_p recovers ordinary PCA.

## Usage

``` r
genpca(
  X,
  A = NULL,
  M = NULL,
  ncomp = NULL,
  method = c("eigen", "spectra", "deflation"),
  constraints_remedy = c("ridge", "error", "clip", "identity"),
  preproc = multivarious::pass(),
  threshold = 1e-06,
  use_cpp = TRUE,
  maxeig = 800,
  warn_approx = TRUE,
  maxit_spectra = 1000,
  tol_spectra = 1e-09,
  verbose = FALSE
)

gmdLA(
  X,
  Q,
  R,
  k = min(n_orig, p_orig),
  n_orig,
  p_orig,
  maxeig = 800,
  tol = 1e-08,
  use_dual = FALSE,
  warn_approx = TRUE,
  verbose = FALSE
)

gmd_deflationR(X, Q, R, k, thr = 1e-06, verbose = FALSE)

# S3 method for class 'genpca'
truncate(x, ncomp)

# S3 method for class 'genpca'
reconstruct(
  x,
  comp = 1:multivarious::ncomp(x),
  rowind = NULL,
  colind = NULL,
  ...
)

# S3 method for class 'genpca'
ncomp(x)
```

## Arguments

- X:

  Numeric matrix n x p.

- A:

  Column constraint: vector (implies diagonal), dense matrix, or sparse
  symmetric p x p PSD matrix. If \`NULL\`, defaults to identity.

- M:

  Row constraint: vector (implies diagonal), dense matrix, or sparse
  symmetric n x n PSD matrix. If \`NULL\`, defaults to identity.

- ncomp:

  Number of components to extract. Defaults to \`min(dim(X))\`. Must be
  positive.

- method:

  Character string specifying the computation method. One of `"eigen"`
  (default, uses `gmdLA`), `"spectra"` (uses matrix-free C++/Spectra
  implementation `gmd_fast_cpp`), or `"deflation"` (uses
  `gmd_deflationR` or `gmd_deflation_cpp`).

- constraints_remedy:

  Character string specifying the remedy for constraints. One of
  `"error"`, `"ridge"`, `"clip"`, or `"identity"`.

- preproc:

  Pre‑processing transformer object from the \*\*multivarious\*\*
  package (default \`multivarious::pass()\`). Use
  \`multivarious::center()\` for centered GPCA. See
  \`?multivarious::prep\` for options.

- threshold:

  Convergence tolerance for the `"deflation"` method's inner loop.
  Default \`1e-6\`.

- use_cpp:

  Logical. If \`TRUE\` (default) and package was compiled with C++
  support, use faster C++ implementation for `method = "deflation"`.
  Fallback to R otherwise. (Ignored for `method = "eigen"` and
  `method = "spectra"`).

- maxeig:

  Upper bound on subspace dimension for eigen/SVD calculations,
  primarily for `method = "eigen"`. If a constraint matrix dimension is
  `<= maxeig` a full eigen decomposition is used. Otherwise only the
  leading `maxeig` eigencomponents are computed via
  [`RSpectra::eigs_sym`](https://rdrr.io/pkg/RSpectra/man/eigs.html), so
  results may be approximate. Default \`800\`.

- warn_approx:

  Logical. If `TRUE` (default) a warning is emitted when an approximate
  eigen decomposition is used because the dimension exceeds `maxeig`.

- maxit_spectra:

  Maximum iterations for the Spectra iterative solver when
  `method = "spectra"`. Default \`1000\`.

- tol_spectra:

  Tolerance for the Spectra iterative solver when `method = "spectra"`.
  Default \`1e-9\`.

- verbose:

  Logical. If \`TRUE\`, print progress messages. Default \`FALSE\`.

## Value

An object of class \`c("genpca", "bi_projector")\` inheriting from
\`multivarious::bi_projector\`, with slots including:

- u,v:

  Left/right singular vectors scaled by the constraint metrics (MU, AV).
  These correspond to components in the original space's geometry. Use
  \`components(fit)\`.

- ou,ov:

  Orthonormal singular vectors in the constraint metric (U, V such that
  UT M U = I, VT AV = I). These are the core mathematical factors.

- sdev:

  Generalised singular values d_k.

- s:

  Scores ( X V or equivalently MU D). Represent projection of rows onto
  components. Use \`scores(fit)\`.

- preproc:

  The \`multivarious\` pre‑processing object used.

- A, M:

  The constraint matrices used (potentially after coercion to sparse
  format).

- propv:

  Proportion of generalized variance explained by each component.

- cumv:

  Cumulative proportion of generalized variance explained.

## Details

\`gmdLA\` caches the eigen decomposition of the constraint matrices by
storing it as an attribute on the matrix. \`compute_sqrtm()\` returns
this modified matrix so callers can reassign it (e.g. \`R \<-
sqrtm_res\$matrix\`) to reuse the cached decomposition in subsequent
calls.

## Method

We compute the rank‑ncomp factors UDVT that minimise \$\$ \\X -
UDV^\top\\\_{M,A}^2 = \mathrm{tr}\\\bigl(M\\
(X-UDV^\top)\\A\\(X-UDV^\top)^\top\bigr) \$\$ subject to UT M U = I, VT
AV = I. (Allen et al., 2014). Three methods are available via the
\`method\` argument:

- `"eigen"` (Default): Uses a one-shot eigen decomposition strategy
  based on `gmdLA`. It explicitly forms and decomposes a \\p \times p\\
  or \\n \times n\\ matrix (depending on `n` vs `p`).

- `"spectra"`: Uses a matrix-free iterative approach via the RcppSpectra
  package to solve the same eigen problem as `"eigen"` but without
  forming the large intermediate matrix. Generally faster and uses less
  memory for large `n` or `p`. Requires C++ compiler and RcppSpectra.

- `"deflation"`: Uses an iterative power/deflation algorithm. Can be
  slower but potentially uses less memory than `"eigen"` for very large
  dense problems where `ncomp` is small.

For pre-computed covariance matrices C = X'MX, see
[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
which performs GPCA directly on C with column constraint R (equivalent
to A).

## References

Allen, G. I., Grosenick, L., & Taylor, J. (2014). \*A Generalized
Least‑Squares Matrix Decomposition.\* Journal of the American
Statistical Association, 109(505), 145‑159. arXiv:1102.3074.

## See also

[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
for GPCA on pre-computed covariance matrices, `truncate.genpca`,
`reconstruct.genpca`, \`multivarious::bi_projector\`,
\`multivarious::project\`, \`multivarious::scores\`,
\`multivarious::components\`, \`multivarious::reconstruct\`.

## Examples

``` r
if (requireNamespace("RSpectra", quietly = TRUE) &&
    requireNamespace("multivarious", quietly = TRUE)) {
  set.seed(123)
  X <- matrix(stats::rnorm(200 * 100), 200, 100)
  rownames(X) <- paste0("R", 1:200)
  colnames(X) <- paste0("C", 1:100)

  # Standard PCA (A=I, M=I, centered) - using default method="eigen"
  gpca_std_eigen <- genpca(X, ncomp = 5, preproc = multivarious::center(), verbose = FALSE)
  
  # Standard PCA using Spectra method (requires C++ build)
  # gpca_std_spectra <- try(genpca(X, ncomp = 5,
  #                              preproc = multivarious::center(),
  #                              method = "spectra", verbose = TRUE))
  # if (!inherits(gpca_std_spectra, "try-error")) {
  #    print(head(gpca_std_spectra$sdev))
  # }

  # Compare singular values with prcomp
  pr_std <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  print("Eigen Method Sdev:")
  print(head(gpca_std_eigen$sdev))
  print("prcomp Sdev:")
  print(head(pr_std$sdev))
  print(paste("Total Var Explained (Eigen):",
              round(sum(gpca_std_eigen$propv) * 100), "%"))

  # Weighted column PCA (diagonal A, no centering)
  col_weights <- stats::runif(100, 0.5, 1.5)
  gpca_weighted <- genpca(X, A = col_weights, ncomp = 3,
                          preproc = multivarious::pass(), verbose = FALSE)
  print("Weighted GPCA Sdev:")
  print(gpca_weighted$sdev)
  print(head(components(gpca_weighted)))
}
#> [1] "Eigen Method Sdev:"
#> [1] 23.38656 23.16950 23.00393 22.56341 22.20862
#> [1] "prcomp Sdev:"
#> [1] 1.657829 1.642442 1.630705 1.599478 1.574327 1.543299
#> [1] "Total Var Explained (Eigen): 13 %"
#> [1] "Weighted GPCA Sdev:"
#> [1] 24.79929 24.23576 23.57776
#> Error in components(gpca_weighted): could not find function "components"
```
