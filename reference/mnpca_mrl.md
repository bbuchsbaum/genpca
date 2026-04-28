# Matrix-Normal PCA via Maximum Regularized Likelihood

Fits a rank-`ncomp` matrix factorization under matrix-normal noise with
sparse row/column precision matrices: \$\$ Y = XW^\top + E,\quad E \sim
\mathcal{MN}(0,\Omega,\Sigma) \$\$ using block coordinate descent:

1.  Alternating least squares updates for `X, W` with fixed precisions
    (`Theta_row = Omega^{-1}`, `Theta_col = Sigma^{-1}`).

2.  Graphical-lasso style precision updates for `Theta_row` and
    `Theta_col` using ADMM, warm starts, and optional block screening.

## Usage

``` r
mnpca_mrl(
  Y,
  ncomp = min(dim(Y)),
  lambda_row = 0.05,
  lambda_col = 0.05,
  max_outer = 25,
  max_inner = 5,
  tol = 1e-04,
  eps_ridge = 1e-08,
  jitter = 1e-06,
  center = TRUE,
  update_precisions = TRUE,
  warm_start = TRUE,
  gl_maxit = 200,
  gl_tol = 1e-04,
  gl_rho = 1,
  penalize_diagonal = FALSE,
  block_screen = TRUE,
  scale_fix = c("trace", "none"),
  sparsify_threshold = 1e-08,
  as_sparse_precision = TRUE,
  verbose = FALSE
)
```

## Arguments

- Y:

  Numeric matrix (`n x p`).

- ncomp:

  Target rank (`r`).

- lambda_row:

  L1 penalty for row precision (`Theta_row`).

- lambda_col:

  L1 penalty for column precision (`Theta_col`).

- max_outer:

  Maximum number of outer BCD iterations.

- max_inner:

  Maximum ALS steps per outer iteration.

- tol:

  Relative tolerance used for ALS and objective convergence checks.

- eps_ridge:

  Ridge added to small `r x r` normal-equation systems.

- jitter:

  Small diagonal jitter used in covariance/precision updates.

- center:

  Logical; center columns of `Y` before fitting.

- update_precisions:

  Logical; if `FALSE`, keeps identity precisions and runs weighted ALS
  only.

- warm_start:

  Logical; warm start precision updates from previous iterate.

- gl_maxit:

  Maximum ADMM iterations per graphical-lasso subproblem.

- gl_tol:

  ADMM convergence tolerance for graphical-lasso subproblems.

- gl_rho:

  ADMM augmented Lagrangian parameter.

- penalize_diagonal:

  Logical; whether to penalize diagonal precision entries in L1 term.
  Default `FALSE`.

- block_screen:

  Logical; use thresholded connected components to solve precision
  subproblems blockwise.

- scale_fix:

  One of `"trace"` or `"none"`. `"trace"` normalizes each precision to
  mean diagonal 1 after updates.

- sparsify_threshold:

  Off-diagonal magnitude threshold used to set tiny precision entries to
  zero after each graphical-lasso solve.

- as_sparse_precision:

  Logical; store precisions as sparse matrices when many entries are
  zero.

- verbose:

  Logical; print iteration diagnostics.

## Value

An object of class `"mnpca_mrl"` with components:

- X, W:

  Estimated low-rank factors (`n x r`, `p x r`).

- Theta_row, Theta_col:

  Estimated row/column precision matrices.

- fitted:

  Reconstructed matrix on input scale.

- fitted_centered:

  Reconstructed centered matrix used in optimization.

- residual_centered:

  Centered residual matrix.

- objective_path:

  Objective values across outer iterations.

- iterations:

  Number of outer iterations used.

- converged:

  Logical convergence flag for outer loop.

- center:

  Column centering vector (or `NULL`).

- call:

  Matched call.

## Details

The covariance updates use low-rank correction identities and avoid
explicit construction of `E = Y - XW^T`.

This implementation follows the maximum regularized likelihood (MRL)
formulation of MN-PCA, combining low-rank factor updates with sparse
precision estimation in row and column spaces.

The main optimization target is: \$\$ \frac12\mathrm{tr}\left(\Theta_c
(Y-XW^\top)^\top \Theta_r (Y-XW^\top)\right)
-\frac{p}{2}\log\|\Theta_r\|-\frac{n}{2}\log\|\Theta_c\| +
n\lambda_r\\\Theta_r\\\_1 + p\lambda_c\\\Theta_c\\\_1 \$\$ with
\\\Theta_r \succ 0, \Theta_c \succ 0\\. When
`update_precisions = FALSE`, the method reduces to weighted low-rank
approximation with fixed identity precisions.

## References

Zhang, C., Gai, K., & Zhang, S. (2024). *Matrix normal PCA for
interpretable dimension reduction and graphical noise modeling*. Pattern
Recognition, 154, 110591.
[doi:10.1016/j.patcog.2024.110591](https://doi.org/10.1016/j.patcog.2024.110591)

Friedman, J., Hastie, T., & Tibshirani, R. (2008). *Sparse inverse
covariance estimation with the graphical lasso*. Biostatistics, 9(3),
432-441.

## Examples

``` r
set.seed(123)
n <- 20; p <- 12; r <- 3
Y <- matrix(rnorm(n * p), n, p)
fit <- mnpca_mrl(
  Y,
  ncomp = r,
  lambda_row = 0.08,
  lambda_col = 0.08,
  max_outer = 6,
  max_inner = 6,
  verbose = FALSE
)
dim(fit$X)
#> [1] 20  3
dim(fit$W)
#> [1] 12  3
length(fit$objective_path)
#> [1] 6
```
