# Sparse and Functional Principal Components Analysis (SFPCA) with Spatial Coordinates

Performs Sparse and Functional PCA on a data matrix, allowing for both
sparsity and smoothness in the estimated principal components. Penalty
parameters left \`NULL\` are selected automatically (see Details). The
spatial smoothness penalty is constructed based on provided spatial
coordinates.

## Usage

``` r
second_diff_matrix(n)

sfpca(
  X,
  K,
  spat_cds,
  lambda_u = NULL,
  lambda_v = NULL,
  alpha_u = NULL,
  alpha_v = NULL,
  Omega_u = NULL,
  penalty_u = "l1",
  penalty_v = "l1",
  nlambda = 10,
  lambda_min_ratio = 0.01,
  knn = min(6, ncol(X) - 1),
  max_iter = 100,
  tol = 1e-06,
  verbose = FALSE,
  uthresh = NULL,
  vthresh = NULL
)

construct_spatial_penalty(spat_cds, method = "distance", k = 6L)
```

## Arguments

- X:

  A numeric data matrix of dimensions n (observations/time points) by p
  (variables/space).

- K:

  The number of principal components to estimate.

- spat_cds:

  A matrix of spatial coordinates for each column of X (variables). Each
  row corresponds to a spatial dimension (e.g., x, y, z), and each
  column corresponds to a variable.

- lambda_u:

  Sparsity penalty parameter for u. If NULL, selected per component by
  BIC along a regularization path (see Details).

- lambda_v:

  Sparsity penalty parameter for v. If NULL, selected per component by
  BIC along a regularization path (see Details).

- alpha_u:

  Smoothness penalty parameter for u. If NULL, defaults to \`1 /
  lambda_max(Omega_u)\` (see Details).

- alpha_v:

  Smoothness penalty parameter for v. If NULL, defaults to \`1 /
  lambda_max(Omega_v)\` (see Details).

- Omega_u:

  A positive semi-definite matrix for smoothness penalty on u. If NULL,
  defaults to second differences penalty (sparse matrix).

- penalty_u:

  The penalty function for u. Either "l1" (lasso) or "scad".

- penalty_v:

  The penalty function for v. Either "l1" (lasso) or "scad".

- nlambda:

  Number of values on the regularization path used for BIC selection of
  \`lambda_u\`/\`lambda_v\` when they are NULL.

- lambda_min_ratio:

  Smallest path value as a fraction of the closed-form \`lambda_max\`,
  on a log-spaced grid.

- knn:

  Number of nearest neighbours for constructing \`Omega_v\`.

- max_iter:

  Maximum number of iterations for the alternating optimization.

- tol:

  Tolerance for convergence of the rank-1 objective.

- verbose:

  Logical; if TRUE, prints progress messages.

- uthresh:

  Deprecated and ignored; \`lambda_u\` is now selected by BIC.

- vthresh:

  Deprecated and ignored; \`lambda_v\` is now selected by BIC.

## Value

An object of class \`c("sfpca", "bi_projector")\` from the multivarious
framework. Use \`multivarious::scores()\` for the sample scores (\\U
D\\), \`multivarious::components()\` for the sparse loadings \\V\\,
\`multivarious::sdev()\` for the singular values, and
\`multivarious::reconstruct()\` for the rank-\`K\` approximation. The
selected penalty parameters are stored as \`lambda_u\`, \`lambda_v\`,
\`alpha_u\`, and \`alpha_v\`. For backward compatibility the pre-0.1
list fields \`\$d\` (singular values) and \`\$u\` (left factors) remain
readable but emit a deprecation warning; use \`sdev()\` and
\`scores()\`/\`\$ou\` instead.

## Details

Each rank-1 problem is solved by alternating solves of the penalized
quadratic subproblems (via C++ coordinate descent) followed by rescaling
onto the smoothness-metric ball, in the constraint form of Allen &
Weylandt (2019). For the convex \`"l1"\` penalty with subproblems solved
to tolerance (the internal \`exact_inner = TRUE\` path, used by the
monotonicity test) the objective is monotonically non-decreasing; the
default inexact path tightens the inner tolerance to a floor before it
may declare convergence, reproducing the same terminal iterates but
without an every-iteration monotonicity guarantee (it may also stop at
\`max_iter\`).

When \`lambda_u\` or \`lambda_v\` is \`NULL\` it is selected per
component by a BIC-style criterion along a regularization path. For the
convex \`"l1"\` penalty \`lambda_max = max(abs(b))\` is, in closed form,
the smallest value whose subproblem solution is exactly zero (at \`x =
0\` the \`S x\` term vanishes, so the KKT condition \`\|b_j\| \<=
lambda\` does not depend on \`S\`); where \`b\` is the matrix-vector
product with the other factor fixed at the SVD initializer. For the
non-convex \`"scad"\` penalty the same value anchors the path but is not
a global-optimality threshold. \`nlambda\` values are laid log-spaced
down to \`lambda_min_ratio \* lambda_max\`, coordinate descent is
warm-started along the path, and the value minimizing \`log(RSS / (n
p)) + df \* log(n p) / (n p)\` is chosen, with \`df\` the support size
of the solution and \`RSS\` the one-sided rank-1 residual sum of squares
with the opposite factor held fixed (a selection heuristic, not the BIC
of the fully alternated rank-1 model). The all-zero solution (at
\`lambda_max\`) is a legitimate candidate: if no rank-1 structure
justifies its degrees of freedom, the component is returned as exactly
zero with \`d = 0\`.

When \`alpha_u\` or \`alpha_v\` is \`NULL\` it defaults to \`1 /
lambda_max(Omega)\`, so the roughest direction of the smoothness penalty
is weighted exactly as strongly as the identity term. This makes the
default invariant to the scaling of \`Omega\` and bounds the condition
number of every subproblem system \`I + alpha \* Omega\` by 2.

## References

Allen, G. I., & Weylandt, M. (2019). Sparse and functional principal
components analysis. In *2019 IEEE Data Science Workshop (DSW)* (pp.
11-16).

## See also

\[genpca()\] for the shared multivarious verbs;
[`multivarious::bi_projector`](https://bbuchsbaum.github.io/multivarious/reference/bi_projector.html).

## Examples

``` r
library(Matrix)
set.seed(123)
# Smooth temporal factor, sparse spatial factor
n <- 100  # Number of time points
p <- 50   # Number of spatial locations
u <- sin(seq(0, 2 * pi, length.out = n))
v <- c(rnorm(10), rep(0, p - 10))
X <- 8 * tcrossprod(u / sqrt(sum(u^2)), v / sqrt(sum(v^2))) +
  matrix(rnorm(n * p, sd = 0.2), n, p)
spat_cds <- matrix(runif(p * 3), nrow = 3, ncol = p)  # 3D coordinates
result <- sfpca(X, K = 1, spat_cds = spat_cds)
multivarious::sdev(result)                  # captured covariance (BIC-tuned)
#> [1] 7.604738
sum(multivarious::components(result) != 0)  # sparse spatial loading
#> [1] 7
```
