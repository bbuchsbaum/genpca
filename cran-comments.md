# genpca 0.2.2 check notes (candidate)

## Reported failure and fix

CRAN's Fedora R-devel checks of 0.2.1 report `starting vector near the null
space` in `test-gplssvd-op-large.R:99`. irlba 2.4.1 removed the `mult` argument;
the old implementation supplied a zero placeholder matrix and a callback
through that argument. The replacement supplies a matrix-free S4 operator
with forward and adjoint multiplication methods using irlba's documented
interface. No tests were skipped, seeds changed, or tolerances loosened.

Outside the existing small-dense fallback, irlba's own tiny-matrix fallback
would require materialization when the smaller dimension is below six. Such
shapes now receive an explicit error directing callers to the eigencore
backend, rather than returning zero singular values or silently densifying.

## Local reproduction and regression results

On Debian 13, R 4.5.0, GCC 14.2.0, Matrix 1.7.3, with all suggested test
dependencies installed:

- Unmodified 0.2.1 + source-installed irlba 2.4.1: 1,582 passing expectations,
  one null-space error at the reported test, one existing empty-test skip.
- Unmodified 0.2.1 + irlba 2.3.7: 1,588 passing expectations, no failures.
- Patched 0.2.2 + irlba 2.4.1: 1,620 passing expectations, no failures/errors.
- Patched 0.2.2 + irlba 2.3.7: 1,620 passing expectations, no failures/errors.

All four runs report the same 73 warning messages (66 associated with tests
and seven top-level deprecation warnings), and one existing empty-test skip.
Warnings cover legacy testthat deprecations, intentional degenerate-rank
examples, and metric-repair diagnostics. The fix adds no warning messages.

The full source build, including regenerated vignettes, passed locally.
Full R CMD check and hosted validation are in progress; these local test
results alone are not a claim of a clean CRAN check.

## Hosted validation and environment limits

Draft pull request with current check links and retained diagnostics:
<https://github.com/bbuchsbaum/genpca/pull/52>.

The workflow source-installs and asserts irlba 2.4.1 (to avoid stale 2.3.7
binary caches), records compiler/R/BLAS and dependency versions, retains
complete warning details and source hashes, and checks:

- R-hub gcc16 / Fedora 44 / R-devel.
- R-hub clang23 / Ubuntu / R-devel.
- R-hub Ubuntu R-release, including a full PDF-manual and vignette check.

R-hub's clang23 host is not CRAN's Fedora Clang host. Its OS, R revision,
Clang patch release and BLAS can differ. R-hub has no native four-way
GCC16/Clang23 by R-devel/R-release matrix. Exact runtime differences and
terminal check results must be reviewed before CRAN submission. No CRAN
submission has been performed.
