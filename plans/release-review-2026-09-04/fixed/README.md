# Release review remediation evidence

The six findings in `../../2026-09-04-release-review.md` have been addressed.
The initial review and its logs remain historical evidence.

- Positive diagonal weights use the same support in forward and inverse
  factors throughout GMD, covariance PCA and PLS. All positive diagonal
  weights are retained; general eigendecompositions retain their documented
  numerical-null-space policy.
- C++ and R deflation take `rank_rtol` separately from the iteration `threshold`.
  The common return boundary filters vectors and variance summaries together.
  C++ empty subviews are materialized so zero-component vectors are empty.
- Shared eigencore wrappers raise `genpca_solver_nonconvergence` before callers
  can consume an incomplete result. Existing dense fallbacks handle this
  condition; matrix-free PLS propagates it without an unbounded dense allocation.
- Explicit clip reaches strict clipping even within the PSD validation margin
  and reports any actual change.
- `geigen_cov()` documents and tests the equation projected onto the retained
  range, including a nonzero residual for the unprojected equation and a case
  where commuting metrics change the ordering of components.
- Manual equations use valid markup. Diagnostic files are excluded through
  `.Rbuildignore`. Vignette setup no longer leaves a global title-check option
  changed. Documentation and Rcpp bindings have been regenerated.

Regression tests: `tests/testthat/test-release-review-regressions.R` (101
expectations). The complete local suite passes 1,403 expectations with zero
failures, 74 existing testthat warnings and one empty-test skip. Those testthat
warnings are distinct from R CMD check status.

Tarball: `/tmp/genpca-release-fixed/genpca_0.2.0.tar.gz`

SHA-256: `3a73c06123b83fb7a61e4df81f77d46150b47ce54e049c4d39a558d33a910335`

The full checks run on this exact tarball, with vignettes rebuilt and PDF/HTML
manual checks enabled, under R 4.5.1 / macOS 14.3 arm64. `Makevars` is a
check-only configuration selecting Apple clang 15 and the existing Homebrew
Fortran runtime. No compiler warning suppression flags were added; the user's
compiler settings and installed eigencore 1.3.0 were not changed. CRAN
1.0.3 is isolated in `/tmp/genpca-release-review/cranlib`.

Commands (run each check in a separate scratch directory):

```sh
R_MAKEVARS_USER=/tmp/genpca-release-fixed/Makevars LC_ALL=C LANG=C \
  R CMD build /Users/bbuchsbaum/code/genpca
R_MAKEVARS_USER=/tmp/genpca-release-fixed/Makevars LC_ALL=C LANG=C \
  R CMD check --as-cran /tmp/genpca-release-fixed/genpca_0.2.0.tar.gz
R_MAKEVARS_USER=/tmp/genpca-release-fixed/Makevars \
  R_LIBS=/tmp/genpca-release-review/cranlib LC_ALL=C LANG=C \
  R CMD check --as-cran /tmp/genpca-release-fixed/genpca_0.2.0.tar.gz
```

Final status for both full checks: **0 errors, 0 warnings, 2 notes**. Both
passed all 1,403 test expectations, examples, vignette rebuilds and the PDF
manual. HTML math errors are gone. The two remaining notes are:

1. CRAN incoming feasibility: new submission and maintainer identity.
2. Local HTML Tidy is too old, so structural HTML validation is skipped.

The source tarball was inspected: no Rplots.pdf, _problems directory,
testthat-problems.rds or plans directory is included. Submission notes now
reflect these receipts and explicitly mark Linux/Windows and external builder
validation as outstanding. No commit, push or CRAN submission was made.

For Git publication, copied console logs have trailing whitespace removed.
Original raw check logs remain in the scratch check directories above.
