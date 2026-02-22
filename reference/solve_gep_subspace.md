# Fast sub-space solver for a small block of generalized eigen-pairs

Uses pre-conditioned sub-space iteration on the operator \\S_2^{-1}
S_1\\ (or its inverse) to obtain the \`q\` largest or smallest
generalized eigen-values/vectors of \\S_1 v = \lambda S_2 v\\.

## Usage

``` r
solve_gep_subspace(
  S1,
  S2,
  q = 2,
  which = c("largest", "smallest"),
  max_iter = 100,
  tol = 1e-06,
  V0 = NULL,
  seed = NULL,
  reg_S = 0.001,
  reg_T = 1e-06,
  verbose = FALSE
)
```

## Arguments

- S1, S2:

  Symmetric positive-(semi)definite \`dgCMatrix\` (or dense) matrices of
  the same dimension \\d\times d\\.

- q:

  Number of eigen-pairs required (\`q \<\< d\`).

- which:

  \`"largest"\` or \`"smallest"\`.

- max_iter, tol:

  Stopping rule - iteration stops when \`max(abs(lambda_new -
  lambda_old)/abs(lambda_old)) \< tol\`.

- V0:

  Optional \`d x q\` initial block (will be orthonormalised).

- seed:

  Optional integer seed for reproducible random initialisation.

- reg_S, reg_T:

  Ridge terms added to \`S1\`/\`S2\` and the small \`q x q\` Gram matrix
  to guarantee invertibility.

- verbose:

  Logical - print convergence info.

## Value

A list with components

- values:

  length-\`q\` numeric vector of Ritz eigen-values.

- vectors:

  \`d x q\` matrix, columns are orthonormal eigen-vectors in the
  \*original\* S-inner-product.
