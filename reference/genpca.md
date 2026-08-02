# Generalised Principal Components Analysis (GPCA)

Implements the Generalised Least-Squares Matrix Decomposition of Allen,
Grosenick & Taylor (2014) for data observed in a **row** inner-product
space M and a **column** inner-product space A. Setting M = I_n, A = I_p
recovers ordinary PCA.

## Usage

``` r
genpca(
  X,
  A = NULL,
  M = NULL,
  ncomp = NULL,
  method = c("eigen", "auto", "spectra", "randomized", "deflation"),
  constraints_remedy = c("ridge", "error", "clip", "identity"),
  preproc = multivarious::pass(),
  threshold = 1e-06,
  maxit_deflation = 500L,
  use_cpp = TRUE,
  maxeig = 800,
  warn_approx = TRUE,
  maxit_spectra = 1000,
  tol_spectra = 1e-09,
  oversample = 20L,
  n_power = 1L,
  n_polish = 0L,
  jitter_metric = 1e-10,
  seed_randomized = 1234L,
  tol_polish_randomized = 1e-04,
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric matrix n x p.

- A:

  Column constraint: vector (implies diagonal), dense matrix, or sparse
  symmetric p x p PSD matrix. If `NULL`, defaults to identity.

- M:

  Row constraint: vector (implies diagonal), dense matrix, or sparse
  symmetric n x n PSD matrix. If `NULL`, defaults to identity.

- ncomp:

  Number of components to extract. Defaults to `min(dim(X))`. Must be
  positive.

- method:

  Character string specifying the computation method. One of `"eigen"`
  (default, uses `gmdLA`), `"auto"` (heuristic choice among `"eigen"`,
  `"spectra"`, and `"randomized"`), `"spectra"` (uses matrix-free
  C++/Spectra implementation `gmd_fast_cpp`), `"randomized"`
  (approximate randomized block solver `gmd_randomized`), or
  `"deflation"` (uses `gmd_deflationR` or `gmd_deflation_cpp`).

- constraints_remedy:

  Character string specifying how a supplied `A` or `M` that is not
  symmetric positive (semi)definite is repaired. Default `"ridge"`. One
  of: `"error"` (reject the input with an error), `"ridge"` (Gershgorin
  diagonal shift: add the smallest diagonal loading that restores
  positive definiteness, falling back to
  [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) for
  small dense matrices), `"clip"` (spectral clip to the PSD cone by
  zeroing negative eigenvalues; this densifies the matrix and refuses
  sparse input larger than 2000 rows/cols, where `"ridge"` should be
  used instead), or `"identity"` (replace the matrix with the identity).
  Note that
  [`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
  defaults to `"error"` instead of `"ridge"`, since it expects an
  already-validated covariance matrix; see its documentation for
  details.

- preproc:

  Pre-processing transformer object from the **multivarious** package
  (default
  [`multivarious::pass()`](https://bbuchsbaum.github.io/multivarious/reference/pass.html)).
  Use
  [`multivarious::center()`](https://bbuchsbaum.github.io/multivarious/reference/center.html)
  for centered GPCA. See
  [`?multivarious::prep`](https://bbuchsbaum.github.io/multivarious/reference/prep.html)
  for options.

- threshold:

  Convergence tolerance for the `"deflation"` method's inner loop.
  Default `1e-6`. Cutoffs are relative to the scale of the problem (the
  norm/singular-value floors scale with \\\sqrt{\mathrm{tr}(X'MXA)}\\),
  so results are invariant to rescaling `X`. The convergence check is on
  a squared step difference, so the resulting singular-vector accuracy
  scales like \\\sqrt{\code{threshold}}\\, not `threshold` itself.

- maxit_deflation:

  Maximum iterations per component for the `"deflation"` method. Default
  `500`.

- use_cpp:

  Logical. If `TRUE` (default) and package was compiled with C++
  support, use faster C++ implementation for `method = "deflation"`.
  Fallback to R otherwise. (Ignored for `method = "eigen"` and
  `method = "spectra"`).

- maxeig:

  Upper bound on subspace dimension for eigen/SVD calculations,
  primarily for `method = "eigen"`. If a constraint matrix dimension is
  `<= maxeig` a full eigen decomposition is used. Otherwise only the
  leading `maxeig` eigencomponents are computed via
  [`RSpectra::eigs_sym`](https://rdrr.io/pkg/RSpectra/man/eigs.html), so
  results may be approximate. Default `800`.

- warn_approx:

  Logical. If `TRUE` (default) a warning is emitted when an approximate
  eigen decomposition is used because the dimension exceeds `maxeig`.

- maxit_spectra:

  Maximum iterations for the Spectra iterative solver when
  `method = "spectra"`. Default `1000`.

- tol_spectra:

  Tolerance for the Spectra iterative solver when `method = "spectra"`.
  Default `1e-9`.

- oversample:

  Oversampling for `method = "randomized"` (sketch size =
  `ncomp + oversample`). Default `20`.

- n_power:

  Number of power iterations for `method = "randomized"`. Default `1`.

- n_polish:

  Number of optional block-polish iterations for
  `method = "randomized"`. Default `0`.

- jitter_metric:

  Small jitter used in metric orthonormalization for
  `method = "randomized"`. Default `1e-10`.

- seed_randomized:

  Optional seed for `method = "randomized"`. Default `1234`. This fully
  determines the randomized backend's random stream: the C++ kernel
  seeds its own generator from this value rather than from R's
  [`set.seed()`](https://rdrr.io/r/base/Random.html)/`.Random.seed`, and
  calling `genpca()` with `method = "randomized"` does not alter the
  caller's `.Random.seed`. To reproduce a randomized fit, fix
  `seed_randomized`, not the R seed.

- tol_polish_randomized:

  Relative tolerance used for early stopping of polish iterations in
  `method = "randomized"`. Set `0` to disable early stop. Default
  `1e-4`.

- verbose:

  Logical. If `TRUE`, print progress messages. Default `FALSE`.

## Value

An object of class `c("genpca", "bi_projector")` inheriting from
[`multivarious::bi_projector`](https://bbuchsbaum.github.io/multivarious/reference/bi_projector.html),
with slots including:

- u,v:

  Left/right singular vectors scaled by the constraint metrics (MU, AV).
  These correspond to components in the original space's geometry. Use
  `components(fit)`.

- ou,ov:

  Orthonormal singular vectors in the constraint metric (U, V such that
  UT M U = I, VT AV = I). These are the core mathematical factors.

- sdev:

  Generalised singular values d_k. Note these are singular values of the
  metric-whitened data matrix, not standard deviations: with identity
  metrics and centering, `sdev = prcomp(X)$sdev * sqrt(nrow(X) - 1)`.

- s:

  Scores: the generalised principal components
  `z_k = X A ov_k = ou_k d_k` (Allen et al. 2014, Section 2.4).
  Identical to `project(fit, X)` on the training data. Use
  `scores(fit)`.

- preproc:

  The `multivarious` pre-processing object used.

- A, M:

  The constraint matrices used (potentially after coercion to sparse
  format).

- propv:

  Proportion of generalized variance explained by each component.

- cumv:

  Cumulative proportion of generalized variance explained.

## Method

We compute the rank-ncomp factors UDVT that minimise \$\$ \\X -
UDV^\top\\\_{M,A}^2 = \mathrm{tr}\\\bigl(M\\
(X-UDV^\top)\\A\\(X-UDV^\top)^\top\bigr) \$\$ subject to UT M U = I, VT
AV = I. (Allen et al., 2014). Five methods are available via the
`method` argument:

- `"eigen"` (Default): Uses a one-shot eigen decomposition strategy
  based on `gmdLA`. It explicitly forms and decomposes a \\p \times p\\
  or \\n \times n\\ matrix (depending on `n` vs `p`).

- `"auto"`: Chooses among `"eigen"`, `"spectra"`, and `"randomized"`
  using heuristics on shape, rank ratio (`ncomp / min(n,p)`), and
  constraint structure.

- `"spectra"`: Uses a matrix-free iterative approach via the RSpectra
  package to solve the same eigen problem as `"eigen"` but without
  forming the large intermediate matrix. Generally faster and uses less
  memory for large `n` or `p`. Requires C++ compiler and RSpectra.

- `"randomized"`: Uses a randomized block range finder and small
  projected eigendecomposition. This is an approximate low-pass method
  that is often much faster for wide dense matrices with sparse metrics
  when only top components are needed.

- `"deflation"`: Uses an iterative power/deflation algorithm. Can be
  slower but potentially uses less memory than `"eigen"` for very large
  dense problems where `ncomp` is small.

## Backend Guidance

The default is `method = "eigen"`; `"auto"` is opt-in, not the default.

- Use `"eigen"` (the default) when you need a stable reference solution
  on small/medium problems.

- Use `"auto"` to let a heuristic pick among `"eigen"`, `"spectra"`, and
  `"randomized"` based on problem shape and constraint structure.

- Use `"spectra"` for larger matrix-free iterative solves where memory
  pressure is a concern.

- Use `"randomized"` for wide low-rank settings (`p >> n`) with sparse
  metrics when throughput matters most.

- Use `"deflation"` when you only need a few components and can tolerate
  iterative convergence behavior.

For pre-computed covariance matrices C = X'MX, see
[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
which performs GPCA directly on C with column constraint R (equivalent
to A).

## References

Allen, G. I., Grosenick, L., & Taylor, J. (2014). *A Generalized
Least-Squares Matrix Decomposition.* Journal of the American Statistical
Association, 109(505), 145-159. arXiv:1102.3074.

## See also

[`genpca_cov`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
for GPCA on pre-computed covariance matrices,
[`truncate.genpca`](https://bbuchsbaum.github.io/genpca/reference/truncate.genpca.md),
[`reconstruct.genpca`](https://bbuchsbaum.github.io/genpca/reference/reconstruct.genpca.md),
[`multivarious::bi_projector`](https://bbuchsbaum.github.io/multivarious/reference/bi_projector.html),
[`multivarious::project`](https://bbuchsbaum.github.io/multivarious/reference/project.html),
[`multivarious::scores`](https://bbuchsbaum.github.io/multivarious/reference/scores.html),
[`multivarious::components`](https://bbuchsbaum.github.io/multivarious/reference/components.html),
[`multivarious::reconstruct`](https://bbuchsbaum.github.io/multivarious/reference/reconstruct.html).

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
#>            PC1         PC2           PC3
#> C1  0.02654255 -0.04444090  6.316309e-02
#> C2 -0.38534643  0.13459026  2.188201e-05
#> C3  0.06053140 -0.07895713  9.358779e-02
#> C4 -0.04387758  0.01422663 -3.160916e-03
#> C5 -0.06087768 -0.01772540  7.120077e-02
#> C6  0.00294182  0.08799976  1.743206e-01
```
