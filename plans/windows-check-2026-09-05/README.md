# Windows check remediation — 2026-09-05

## Report and reproduction

The supplied Windows Server 2022 / R-devel r90492 / GCC 14.3.0 check reported
1 ERROR and 1 NOTE: 10 failed expectations, 84 warnings, 1 skip, 1393 passes.
The source hash of that remote upload was not supplied. No new Windows run
has been performed for the corrected candidate.

The randomized failure was reproduced on macOS by changing the sketch seed
from 1234 to 1 (and several others). With X = diag(1, 1e6) and a row or column
metric diag(1, 1e-10), the analytic singular values are 10 and 1. The old
solver lost rank or returned inconsistent factors. The sign sketch can have
duplicate columns; C++ standard-library distribution implementations do not
promise identical seed-to-sketch mappings. Regularizing a deficient Gram
then incorrectly treated it as a full basis. See before.log and reproduce.R.

The exact Windows likelihood discrepancy was not reproduced locally at
seed 2 / scale 100, but other seeds at that scale exceeded the same bound.
Larger-scale probes also exposed larger final reevaluation discrepancies.
The diagnostic combined reciprocal rescaling with final refitting and
objective reevaluation, contradicting its promise of zero when no rescale
was requested. Removing the precision-to-covariance inverse round trip alone
did not eliminate the discrepancy and was not adopted as the fix.

## Changes

- Full-dimensional sketches use deterministic complete coordinates; partial
  sketches use Gaussian draws, retaining seeded RNG isolation.
- Both R and C++ orthonormalizers validate the candidate in the original
  metric and correct its residual normalization error. Deficient sketches
  use column-scaled SVD followed by metric normalization on the independent
  span. Gram jitter no longer invents usable rank. The component cutoff is
  applied to decomposition singular values, not to the internal metric Gram.
- `loglik_rescale_delta` measures only the penalty effect of reciprocal
  rescaling. It is exactly zero for `scale_fix = "none"`.
  `loglik_refit_delta` retains the remaining final-refit/reevaluation change;
  it is not clamped or concealed. `loglik`, metric learning, convergence, and
  the objective-at-returned-metrics checks remain unchanged.
- Tests compare analytic spectra, both metric orthogonality identities,
  reconstructions, multiple seeds, R/C++, all dense/sparse metric dispatches,
  rank-deficient blocks, partial sketches, and a direct likelihood oracle.
  The old large-scale assertion now tests the correctly separated quantities;
  the monotonicity and positive-definiteness assertions were retained.
- Updated API help, NEWS, the metric vignette, and submission comments.
  The incoming spelling NOTE concerns citation author names, already present
  in inst/WORDLIST; their spelling is explained in cran-comments.md.

## Verification

- Full local test suite: 1588 passed expectations, zero failures/errors,
  74 emitted existing testthat warnings (including file-level deprecations),
  and one existing empty-test skip. These are distinct from check status.
- Fresh exact-tarball `R CMD check --as-cran` with eigencore 1.0.3 and 1.3.0:
  both 0 errors, 0 warnings, 2 notes. Tests, examples, vignette rebuilds,
  PDF manual, and HTML manual generation passed. Notes: new submission;
  local HTML Tidy is too old for structural validation.
- R 4.5.1, macOS arm64, Apple clang 15 via the same scratch Makevars as the
  previous release check. No installed dependency or user setting changed.
- Checked tarball sources match current R, C++, Rd, test, and vignette files.
  `git diff --check` passed.

Tarball: `/tmp/genpca-windows-fix/genpca_0.2.0.tar.gz`

SHA-256: `34a8b4123b5ab5eb7c89d5e6f6d8489a2403ec9606498c79531bdfb94eaeac57`

Commands (checks ran in separate scratch directories):

```sh
LC_ALL=C LANG=C R_MAKEVARS_USER=/tmp/genpca-release-fixed/Makevars \
  R CMD build /Users/bbuchsbaum/code/genpca
LC_ALL=C LANG=C R_MAKEVARS_USER=/tmp/genpca-release-fixed/Makevars \
  R CMD check --as-cran /tmp/genpca-windows-fix/genpca_0.2.0.tar.gz
LC_ALL=C LANG=C R_MAKEVARS_USER=/tmp/genpca-release-fixed/Makevars \
  R_LIBS=/tmp/genpca-release-review/cranlib \
  R CMD check --as-cran /tmp/genpca-windows-fix/genpca_0.2.0.tar.gz
```

Logs copied here have trailing whitespace stripped for Git; raw outputs and
check installations remain under `/tmp/genpca-windows-fix`. No commit, push,
external builder upload, or CRAN submission was performed in this step.
