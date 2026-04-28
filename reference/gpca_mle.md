# Experimental ML estimation of GPCA metrics

Alternates between GPCA factor estimation and maximum-likelihood updates
of the row/column metric matrices (M, A) under a Gaussian matrix-normal
error model. Each iteration:

1.  Fit GPCA with current `A`, `M`.

2.  Reconstruct \\\hat X\\; compute residuals \\E = X - \hat X\\.

3.  Update covariance estimates \\\Sigma_r = E A E^T / p\\ and
    \\\Sigma_c = E^T M E / n\\, apply ridge `lambda`, fix scale, set
    `M <- solve(Sigma_r)`, `A <- solve(Sigma_c)`, then project back to
    SPD via `ensure_spd`.

## Usage

``` r
gpca_mle(
  X,
  ncomp = min(dim(X)),
  max_iter = 20,
  lambda = 0.001,
  scale_fix = c("trace", "det", "none"),
  tol = 1e-04,
  method = "eigen",
  constraints_remedy = "ridge",
  preproc = multivarious::pass(),
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  Numeric matrix (n x p).

- ncomp:

  Rank to extract at each GPCA step.

- max_iter:

  Maximum outer alternations (default 20).

- lambda:

  Ridge term added to \\\Sigma_r\\ and \\\Sigma_c\\ before inversion to
  keep them SPD (default 1e-3).

- scale_fix:

  How to resolve the `c * Sigma_r, Sigma_c / c` indeterminacy. One of
  `"trace"` (default, mean diagonal = 1), `"det"` (determinant = 1), or
  `"none"`.

- tol:

  Relative tolerance on successive log-likelihood change (default 1e-4)
  for early stopping.

- method:

  GPCA method passed to `genpca` (defaults to "eigen").

- constraints_remedy:

  Passed to `genpca`; defaults to "ridge".

- preproc:

  Pre-processing transformer; defaults to
  [`multivarious::pass()`](https://bbuchsbaum.github.io/multivarious/reference/pass.html).

- verbose:

  Logical; if TRUE, prints iteration diagnostics.

- ...:

  Additional arguments forwarded to `genpca`.

## Value

A list with elements `fit` (the final `genpca` result), `A`, `M`
(learned SPD metrics), `loglik` (final log-likelihood), and
`loglik_path`.

## Details

This is an experimental convenience wrapper; for stability it enforces
SPD at every step, applies ridge shrinkage, and fixes the scale
indeterminacy via the mean diagonal (`scale_fix = "trace"`).

The matrix-normal likelihood is flat along \\c\\\Sigma_r, \Sigma_c/c\\;
a scale constraint (trace or determinant) is therefore imposed before
inversion. The algorithm stops when the relative change in
log-likelihood is below `tol` or `max_iter` is reached. Increase
`lambda` or reduce `ncomp` if iterations become unstable.

## References

Dutilleul, P. (1999). \*The MLE algorithm for the matrix normal
distribution\*. Journal of Statistical Computation and Simulation,
64(2), 105-123.

## Examples

``` r
if (requireNamespace("multivarious", quietly = TRUE)) {
  set.seed(123)
  X <- matrix(rnorm(40), 8, 5)
  res <- gpca_mle(X, ncomp = 2, max_iter = 5, lambda = 1e-3,
                  scale_fix = "trace", verbose = FALSE)
  # Learned metrics are SPD and match dimensions
  dim(res$A); dim(res$M)
  res$loglik_path
}
#> Warning: `reverse_transform()` was deprecated in multivarious 0.3.0.
#> ℹ Please use `inverse_transform()` instead.
#> ℹ reverse_transform() is deprecated. Use inverse_transform() for a more
#>   standard interface.
#> ℹ The deprecated feature was likely used in the genpca package.
#>   Please report the issue at <https://github.com/bbuchsbaum/genpca/issues>.
#> 'as(<ddenseMatrix>, "dgeMatrix")' is deprecated.
#> Use 'as(as(., "generalMatrix"), "unpackedMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
#> [1] 105.03314  84.68565  90.89784  84.59021  90.72233
```
