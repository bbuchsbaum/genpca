# Canonical Generalized PLS (alias)

Convenience alias for
[`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md);
computes canonical generalized PLS (PLS-SVD/GPLSSVD). See
[`?genpls`](https://bbuchsbaum.github.io/genpca/reference/genpls.md) for
full documentation.

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

  Numeric or Matrix, n x q. Must have same n as `X`.

- Ax:

  Column metric for X (W_X): vector/diagonal/matrix; `NULL` means
  identity.

- Ay:

  Column metric for Y (W_Y): vector/diagonal/matrix; `NULL` means
  identity.

- Mx:

  Row metric for X (M_X): vector/diagonal/matrix; `NULL` means identity.

- My:

  Row metric for Y (M_Y): vector/diagonal/matrix; `NULL` means identity.

- ncomp:

  Number of components to extract (rank-k). Default 2.

- preproc_x, preproc_y:

  Optional `multivarious` preprocessors (e.g., `center()`). Defaults to
  [`multivarious::pass()`](https://bbuchsbaum.github.io/multivarious/reference/pass.html)
  (no-op).

- svd_backend:

  Character, one of `"RSpectra"` (default) or `"irlba"` for the
  iterative SVD. This choice only matters for larger problems: whenever
  both `X` and `Y` have at most 64 columns after preprocessing, the
  operator materializes `S` densely and computes a direct
  [`svd()`](https://rdrr.io/r/base/svd.html), ignoring `svd_backend`
  entirely (see
  [`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md)).

- svd_opts:

  List of options passed to the SVD backend, e.g., `tol`, `maxitr`.

- verbose:

  Logical; print brief progress messages.

## Value

An object of class `c("genpls", "cross_projector", "projector")` with
the same structure as
[`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md)
returns (X-/Y-weights `vx`/`vy`, singular values `d`, generalized
weights `p`/`q`, scores `fi`/`fj`, latent variables `lx`/`ly`, `ncomp`,
and `backend`); see
[`?genpls`](https://bbuchsbaum.github.io/genpca/reference/genpls.md) for
the definition of each slot.

## References

Beaton, D. (2020). Generalized eigen, singular value, and partial least
squares decompositions: The GSVD package. (Eqs. 10-14).
arXiv:2010.14734.

## See also

[`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(60 * 5), 60, 5)
Y <- matrix(rnorm(60 * 4), 60, 4)
fit <- genplsc(X, Y, ncomp = 2,
               preproc_x = multivarious::center(),
               preproc_y = multivarious::center())
fit$d
#> [1] 24.43306 23.88050
```
