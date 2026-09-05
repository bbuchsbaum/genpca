# Utilities for constraints

Helpers to validate, symmetrize and (only when asked) repair constraint
matrices. Two relative tolerances are used throughout the package:
`metric_rtol` (default `sqrt(.Machine$double.eps)`) decides whether a
metric is positive (semi)definite and which of its eigenvalues count as
zero, and `rank_rtol` (see
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md))
decides which components are kept. Both are relative to the scale of the
matrix, so every decision is invariant to rescaling the input.
