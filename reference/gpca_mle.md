# Experimental penalized-ML estimation of GPCA metrics

Alternates between GPCA factor estimation and penalized
maximum-likelihood updates of the row/column metric matrices (M, A)
under a Gaussian matrix-normal error model with a low-rank mean. Each
iteration performs three exact block minimizations of a single penalized
objective:

1.  Fit GPCA with current `A`, `M`: the GMD theorem makes this the best
    rank-`ncomp` fit in the (M, A) norm, so it exactly minimizes the
    residual term.

2.  Update \\\Sigma_r = E A E^T / p + \lambda I\\ (the exact block
    minimizer under the ridge penalty), set `M = solve(Sigma_r)`.

3.  Update \\\Sigma_c = E^T M E / n + \lambda I\\ using the *updated*
    `M` (sequential flip-flop, Dutilleul 1999), set
    `A = solve(Sigma_c)`.

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

  Ridge penalty weight (default 1e-3). Part of the objective (MAP
  interpretation), not just a numerical safeguard: it shrinks both
  covariances toward a multiple of the identity and pins the row/column
  scale split during iteration. Must be non-negative; with `lambda = 0`
  the objective loses strict convexity in the scale direction and
  covariances may become singular.

- scale_fix:

  How to canonicalize the `c * Sigma_r, Sigma_c / c` indeterminacy at
  exit. One of `"trace"` (default: row covariance scaled to mean
  diagonal 1), `"det"` (row covariance scaled to determinant 1), or
  `"none"`. Applied as a joint reciprocal rescale, so the fitted
  covariance \\\Sigma_r \otimes \Sigma_c\\ and the likelihood are
  unchanged.

- tol:

  Relative tolerance on successive penalized log-likelihood change
  (default 1e-4) for early stopping.

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

A list with elements `fit` (a `genpca` fit computed with the returned
canonicalized metrics), `A`, `M` (learned SPD metrics), `loglik` (final
penalized log-likelihood), and `loglik_path` (the penalized
log-likelihood after each outer iteration; monotone non-decreasing up to
numerical noise, since every block update exactly minimizes the shared
penalized objective). Values omit additive constants and include the
`lambda` penalty, so they are comparable across iterations and across
runs with the same `lambda`, but not across different `lambda` values.

## Details

The objective is the matrix-normal log-likelihood with a low-rank mean
and an inverse-Wishart-style ridge penalty
\\\lambda\\(p\\\mathrm{tr}\\\Sigma_r^{-1} +
n\\\mathrm{tr}\\\Sigma_c^{-1})\\ (a MAP estimate). Because every block
update is an exact minimizer of this one objective, `loglik_path` is
monotone non-decreasing up to numerical noise. The penalty also resolves
the \\c\\\Sigma_r, \Sigma_c/c\\ scale indeterminacy during iteration;
`scale_fix` is applied once at exit as a *joint* reciprocal rescale (row
covariance normalized, factor absorbed into the column covariance),
which leaves the likelihood unchanged. The algorithm stops when the
relative change in the penalized log-likelihood falls below `tol` or
`max_iter` is reached. Increase `lambda` or reduce `ncomp` if iterations
become unstable.

## References

Dutilleul, P. (1999). *The MLE algorithm for the matrix normal
distribution*. Journal of Statistical Computation and Simulation, 64(2),
105-123.

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
#> [1] 120.8615 122.9804 125.0796 127.1557 129.2007
```
