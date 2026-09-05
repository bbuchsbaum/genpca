# Repair a metric matrix explicitly

Returns a positive (semi)definite version of `A` together with a
diagnostic report of what was done. This is the explicit counterpart of
the `constraints_remedy` argument of
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md):
nothing in the package repairs a metric silently, and this function lets
you inspect the repair before using the result.

## Usage

``` r
repair_metric(
  A,
  method = c("ridge", "clip", "identity"),
  rtol = .metric_rtol_default(),
  name = "A",
  diag_maxn = 2000L
)
```

## Arguments

- A:

  A square symmetric matrix (base matrix or `Matrix`). Asymmetry beyond
  roundoff is an error (see the relative asymmetry test in
  [`symmetrize_or_stop()`](https://bbuchsbaum.github.io/genpca/reference/symmetrize_or_stop.md));
  it is not something a PSD repair should hide.

- method:

  `"ridge"` adds a diagonal loading that makes the matrix positive
  definite (Gershgorin-based shift, falling back to
  [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) for
  small dense matrices); `"clip"` projects onto the PSD cone by zeroing
  negative eigenvalues (dense eigendecomposition; refuses large sparse
  input); `"identity"` replaces an indefinite matrix by the identity
  (the report still describes the input).

- rtol:

  Relative tolerance: eigenvalues above `-rtol * scale(A)` count as
  non-negative for `"ridge"` and `"identity"`, which then return the
  matrix unchanged. `"clip"` always removes negative eigenvalues,
  regardless of this tolerance (up to reconstruction roundoff). Default
  `sqrt(.Machine$double.eps)`.

- name:

  Label used in messages.

- diag_maxn:

  Largest dimension for which the report computes the full spectrum
  (minimum eigenvalue, rank, condition number); above it an iterative
  minimum eigenvalue estimate and a Gershgorin bound are reported, with
  rank and condition number unavailable.

## Value

The repaired matrix (a `Matrix`), with attribute `"repair_report"` of
class `"metric_repair_report"`: a list with `name`, `method`, `changed`,
`n`, `min_eigenvalue_before`, `min_eigenvalue_after`,
`gershgorin_bound_before`, `shift` (diagonal loading added by `"ridge"`,
`NA` for `"clip"`), `rank`, `condition_number` and `rtol`.

## See also

[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
(argument `constraints_remedy`)

## Examples

``` r
A <- matrix(c(1, 2, 2, 1), 2)           # eigenvalues 3 and -1
B <- repair_metric(A, method = "ridge")
attr(B, "repair_report")
#> Metric repair report for A (2 x 2)
#>   method:                ridge
#>   changed:               TRUE
#>   min eigenvalue before: -1
#>   min eigenvalue after:  0.1
#>   Gershgorin bound:      -1
#>   diagonal shift:        1.1
#>   rank:                  2
#>   condition number:      41
#>   relative tolerance:    1.49e-08
C <- repair_metric(A, method = "clip")
eigen(as.matrix(C))$values
#> [1] 3 0
```
