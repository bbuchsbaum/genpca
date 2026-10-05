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

## Full package checks

Package-content commit: `1a90a8c6bb55cd9c9923e4571f78fa21442b68a7`.
The final metadata-only update to this file does not alter the built package
(`cran-comments.md` is excluded by `.Rbuildignore`).

R-hub run on 2026-10-05:
<https://github.com/bbuchsbaum/genpca/actions/runs/37342894925>.

- Fedora 44, R-devel r90638, GCC/GFortran 16.2.1: Status OK.
- Ubuntu 26.04.1, R-devel r90631, Clang/Flang 23.1.3: Status OK.
- Ubuntu 24.04.5, R 4.6.1, GCC/GFortran 13.3.0: Status OK.
- Additional full R-release `R CMD check --as-cran --no-stop-on-test-error`,
  including PDF manual and regenerated vignettes: Status OK.

All hosted checks have zero errors, check-level warnings and notes. Each
standard check passed 1,620 expectations with one pre-existing empty-test
skip. The workflow reinstalls irlba 2.4.1 from its checksummed CRAN source and
asserts its loaded version after dependency setup, avoiding stale binaries.

Local Debian 13 / R 4.5.0 full `R CMD check --as-cran
--no-stop-on-test-error`: zero errors, zero warnings, one NOTE:
`unable to verify current time`. Incoming feasibility, examples, all tests,
rebuilt vignettes, PDF manual, HTML validation/math rendering and PDF size
checks passed. The exact checked installed artifact also passed all 1,620
expectations with irlba 2.3.7.

## Warning review

Test warning counts are not check-level warnings. Hosted `R CMD check`
reports 83 test warnings; the separately logged test pass with all installed
packages visible reports the original 73. A controlled reproduction using
R CMD check's recommended-package hiding mechanism reproduces exactly 83:
all original 73 warning records are unchanged, plus ten multivarious
regularized-inverse fallback warnings when MASS is unavailable. Hosted
CheckReporter suppresses individual warning records; this attribution is
based on that controlled reproduction and dependency source inspection.
The new operator tests emit no warnings. No diagnostic flags are suppressed.

Clang emits ten compiler warnings from upstream RcppEigen headers
(`-Wunused-but-set-variable`), matching the locations in CRAN's existing
0.2.1 Clang install log. They are retained in the complete install log and
are separate from the test warning counts.

## Environment differences and retained diagnostics

R-hub GCC matches CRAN's current R revision and GCC version, but uses
OpenBLAS 0.3.29 rather than CRAN's reference BLAS. R-hub Clang is Ubuntu,
R-devel r90631, Clang/Flang 23.1.3 and OpenBLAS 0.3.32; CRAN's Clang host is
Fedora, r90638 and Clang/Flang 23.1.2. R-hub's release host uses OpenBLAS
0.3.26. There is no native four-way GCC16/Clang23 by R-devel/R-release matrix.
These checks therefore provide compiler-family and release-R coverage,
not an assertion that every CRAN host has been reproduced exactly.

On the Clang container, igraph's dependency build initially failed because
Flang rejects `-fvisibility=hidden`. A temporary `F_VISIBILITY =` override
is used only while installing dependencies. genpca and the explicitly
reinstalled irlba use the original R-hub compiler settings. The effective
dependency override, original package compiler flags and all versions are
retained in the artifacts. No package compiler warnings are disabled.

Complete check trees, examples/vignettes/manual logs, compiler logs, individual
test warnings/results, environment metadata, timestamps, checked commit and
tarball hashes are attached to the run above. Draft PR and latest checks:
<https://github.com/bbuchsbaum/genpca/pull/52>.

## Checked source archive hashes (SHA-256)

The archives were built from the same package-content commit; generated
vignette outputs and build metadata mean their bytes differ by environment.

- Local full check: `31bfd7a77c0aed807dbceff9a0c5d8f10d4f82e6699bcec832ac90d946dfecda`
- R-hub GCC16: `a49c4ad490deacf263d917e7cc6f0ab457eee16283ed083beda746f5cc856cdd`
- R-hub Clang23: `8b4da7f4eed7d73c1766ce1fe42641f3c8dfc565ea0746511e19cc6a97f98859`
- Full R-release check: `ba2ae46c8fd42ce2ccc3f4409b09961dee7286f2205e02399c69cf1267b59bfd`

No merge or CRAN submission has been performed.
