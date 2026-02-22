# Canonical Generalized PLS (alias)

Convenience alias for \`genpls()\`; computes canonical generalized PLS
(PLS-SVD/GPLSSVD). See \`?genpls\` for full documentation.

## Usage

``` r
genplsc(
  X,
  Y,
  Ax = NULL,
  Ay = NULL,
  Mx = NULL,
  My = NULL,
  ncomp = 2,
  preproc_x = multivarious::pass(),
  preproc_y = multivarious::pass(),
  svd_backend = c("RSpectra", "irlba"),
  svd_opts = list(tol = 1e-07, maxitr = 1000),
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric or Matrix, n x p.

- Y:

  Numeric or Matrix, n x q. Must have same n as \`X\`.

- Ax:

  Column metric for X (W_X): vector/diagonal/matrix; \`NULL\` ⇒
  identity.

- Ay:

  Column metric for Y (W_Y): vector/diagonal/matrix; \`NULL\` ⇒
  identity.

- Mx:

  Row metric for X (M_X): vector/diagonal/matrix; \`NULL\` ⇒ identity.

- My:

  Row metric for Y (M_Y): vector/diagonal/matrix; \`NULL\` ⇒ identity.

- ncomp:

  Number of components to extract (rank-k). Default 2.

- preproc_x, preproc_y:

  Optional \`multivarious\` preprocessors (e.g., \`center()\`). Defaults
  to \`multivarious::pass()\` (no-op).

- svd_backend:

  Character, one of \`"RSpectra"\` (default) or \`"irlba"\` for
  iterative SVD. If neither backend is available, a dense fallback is
  used for small problems by materializing S.

- svd_opts:

  List of options passed to the SVD backend, e.g., \`tol\`, \`maxitr\`.

- verbose:

  Logical; print brief progress messages.

## Value

See \`genpls()\`
