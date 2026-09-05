# Prepare and validate constraint matrices

Coerces `A`/`M` (NULL, weight vector, diagonal, dense or sparse matrix)
to `Matrix` objects and validates them: finite entries, symmetry within
roundoff (see
[`symmetrize_or_stop()`](https://bbuchsbaum.github.io/genpca/reference/symmetrize_or_stop.md)),
and positive semi-definiteness within `tol` relative to the scale of the
matrix. The requested `remedy` is applied to a metric that fails the PSD
check (or has negative eigenvalues within tolerance when explicit
clipping is requested), and every repair emits a warning of class
`genpca_metric_repaired` carrying the
[`repair_metric()`](https://bbuchsbaum.github.io/genpca/reference/repair_metric.md)
report; valid PSD metrics, singular ones included, pass through under
every remedy. Explicit clipping removes even negative eigenvalues within
the validation tolerance. Asymmetric input is an error under every
remedy.

## Usage

``` r
prep_constraints(
  X,
  A,
  M,
  tol = .metric_rtol_default(),
  remedy = c("error", "ridge", "clip", "identity"),
  verbose = FALSE
)
```

## Arguments

- X:

  data matrix (only its dimensions are used)

- A, M:

  column/row constraints

- tol:

  relative PSD tolerance (default `sqrt(.Machine$double.eps)`)

- remedy:

  what to do with an indefinite metric

- verbose:

  emit a message when a metric is replaced by the identity

## Value

list with elements `A` and `M`
