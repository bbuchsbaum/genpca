# Regularised / Generalised Partial Least Squares (RPLS / GPLS)

Implements the algorithm of Allen \*et al.\* (2013) for supervised
dimension-reduction with optional sparsity (\\\ell_1\\) or ridge
(\\\ell_2\\) penalties \*\*and\*\* the generalised extension that
operates in a user-supplied quadratic form \\Q\\.

## Usage

``` r
fit_rpls(
  X,
  Y,
  K = 2,
  lambda = 0.1,
  penalty = c("l1", "ridge"),
  Q = NULL,
  nonneg = FALSE,
  tol = 1e-06,
  maxiter = 200,
  verbose = FALSE
)

rpls(
  X,
  Y,
  K = 2,
  lambda = 0.1,
  penalty = c("l1", "ridge"),
  Q = NULL,
  nonneg = FALSE,
  preproc_x = multivarious::pass(),
  preproc_y = multivarious::pass(),
  tol = 1e-06,
  maxiter = 200,
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  Numeric matrix \\(n \times p)\\ — predictors.

- Y:

  Numeric matrix \\(n \times q)\\ — responses.

- K:

  Integer, number of latent factors to extract. Default \`2\`.

- lambda:

  Scalar or length-`K` numeric vector of penalties.

- penalty:

  Either \`"l1"\` (lasso) or \`"ridge"\`.

- Q:

  Optional positive-(semi)definite \\p \times p\\ matrix inducing
  \*generalised\* PLS. \`NULL\` ⇒ identity.

- nonneg:

  Logical, force non-negative loadings when `penalty = "l1"`. Note: This
  option is currently ignored when `penalty = "ridge"`.

- tol:

  Relative tolerance for the inner iterations convergence check. Default
  \`1e-6\`.

- maxiter:

  Maximum number of inner iterations per component. Default \`200\`.

- verbose:

  Logical; print progress messages during component extraction. Default
  \`FALSE\`.

- preproc_x, preproc_y:

  Optional multivarious preprocessing objects (see
  [`fit_transform`](https://bbuchsbaum.github.io/multivarious/reference/fit_transform.html)).
  By default they pass the data through unchanged using
  [`pass()`](https://testthat.r-lib.org/reference/fail.html).

- ...:

  Further arguments (e.g., custom stopping criteria if implemented) are
  stored in the returned object (they are not used by `fit_rpls`).

## Value

An object of class `c("rpls","cross_projector","projector")` with at
least the elements

- vx:

  \\p \times K\\ matrix of X-loadings.

- vy:

  \\q \times K\\ matrix of Y-loadings.

- ncomp:

  Number of components actually extracted (may be \< K).

- penalty:

  Penalty type used (\`"l1"\` or \`"ridge"\`).

- preproc_x, preproc_y:

  Pre-processing transforms used.

- ...:

  Other parameters like \`lambda\`, \`tol\`, \`maxiter\`, \`nonneg\`,
  \`Q\` indicator, \`verbose\` are also stored.

The object supports [`predict()`](https://rdrr.io/r/stats/predict.html),
`project()`,
[`transfer()`](https://bbuchsbaum.github.io/multivarious/reference/transfer.html),
[`coef()`](https://rdrr.io/r/stats/coef.html) and other multivarious
generics.

## Details

Unlike \`genpls()\` from `genplsr.R`, which handles separate row and
column metrics (\`Mx\`, \`Ax\`, \`My\`, \`Ay\`) with a Gram–Schmidt
orthogonalisation step, \`rpls()\` uses a single metric \`Q\` and the
simpler penalised updates of Allen et al.

## Method

The routine follows Algorithm 1 of Allen \*et al.\* (2013, \*Stat. Anal.
Data Min.\*, 6 : 302–314) — see the paper for details. Briefly, each
component maximises \$\$\max\_{u,v}\\ v^\top Q M u - \lambda \\ P(v)\$\$
with \\Q = I_p\\ for standard RPLS. The alternating updates are: \\u
\leftarrow M^\top Q v / \\M^\top Q v\\\_2\\, then a penalised (possibly
non-negative) regression for \\v\\, normalised in the \\Q\\-norm.

## References

Allen, G. I., Peterson, C., Vannucci, M., & Maletić-Savatić, M. (2013).
\*Regularized Partial Least Squares with an Application to NMR
Spectroscopy.\* \*\*Statistical Analysis and Data Mining, 6(4)\*\*,
302-314. DOI:10.1002/sam.11169.

## Examples

``` r
# Generate sample data
set.seed(123)
n <- 50
p <- 20
q <- 10
X <- matrix(rnorm(n * p), n, p)
Y <- X[, 1:5] %*% matrix(rnorm(5 * q), 5, q) + matrix(rnorm(n * q), n, q)

# Fit regularized PLS with L1 penalty
fit_l1 <- rpls(X, Y, K = 3, lambda = 0.1, penalty = "l1")
print(fit_l1)
#> cross projector:  rpls cross_projector projector 
#> input dim (X):  20 
#> output dim (X):  3 
#> input dim (Y):  10 
#> output dim (Y):  3 

# Fit regularized PLS with ridge penalty
fit_ridge <- rpls(X, Y, K = 3, lambda = 0.1, penalty = "ridge")
print(fit_ridge)
#> cross projector:  rpls cross_projector projector 
#> input dim (X):  20 
#> output dim (X):  3 
#> input dim (Y):  10 
#> output dim (Y):  3 
```
