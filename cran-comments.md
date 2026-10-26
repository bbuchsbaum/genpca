# genpca 0.2.2 CRAN resubmission

## Submission type

This is a corrective resubmission following archival of genpca 0.2.1 on
2026-10-26. It addresses the irlba 2.4.1 compatibility failure reported by
CRAN on 2026-10-05. The correction was merged in PR #52 on 2026-10-05 as
`8f2af8b689da7e986c76d4311bf85df2f0eb0925`. Please consider restoring the
package to CRAN.

## Reported failure and fix

CRAN's Fedora R-devel checks of 0.2.1 reported `starting vector near the null
space` in `test-gplssvd-op-large.R:99`. irlba 2.4.1 removed the `mult` argument;
the old implementation supplied a zero placeholder matrix and a callback
through that argument. The replacement supplies a matrix-free S4 operator
with forward and adjoint multiplication methods using irlba's documented
interface. No tests were skipped, seeds changed, or tolerances loosened.

Outside the existing small-dense fallback, irlba's own tiny-matrix fallback
would require materialization when the smaller dimension is below six. Such
shapes now receive an explicit error directing callers to the eigencore
backend, rather than returning zero singular values or silently densifying.

## Fresh checks (2026-10-26)

Checked source commit: `0421e66723fbf83b9ef6cfd32b9d9fac96239281`.
This differs from merged master only in this file, which is excluded by
`.Rbuildignore`. No package code, documentation, tests, or vignettes changed.

Fresh workflow and complete diagnostics:
<https://github.com/bbuchsbaum/genpca/actions/runs/40326421125>.

| Environment | Result |
| --- | --- |
| Fedora 44, R-devel r90763, GCC/GFortran 16.2.1 | 0 errors, 0 warnings, 0 notes |
| Ubuntu 26.04.1, R-devel r90763, Clang/Flang 23 | 0 errors, 0 warnings, 0 notes |
| Ubuntu 26.04.1, R 4.6.1, GCC/GFortran 15.2.0 | 0 errors, 0 warnings, 0 notes |
| Additional full R 4.6.1 check, including PDF/HTML manual and regenerated vignettes | 0 errors, 0 warnings, 0 notes |
| Local macOS 14.3 arm64, R 4.5.1, Apple clang 15.0.0; full check of local build | 0 errors, 0 warnings, 3 notes |
| Same local environment, full check of exact selected submission archive | 0 errors, 0 warnings, 3 notes |

The hosted standard checks used `--as-cran --no-manual
--no-stop-on-test-error`. The full R-release and local checks used
`--as-cran --no-stop-on-test-error`, without disabling manual or vignette
checks. These environments provide compiler-family and R-version coverage;
they are not exact replicas of every CRAN host.

All environments used CRAN irlba 2.4.1, eigencore 1.0.3 and multivarious
0.3.2. The R-devel and local checks used current CRAN Matrix 1.7-6; the
R-release image supplied recommended Matrix 1.7-5. The hosted workflow
reinstalls irlba from its checksummed CRAN source and verifies its loaded
version. Local updates were installed in an isolated library; the
maintainer's existing library was unchanged.

Each check passed all 1,620 expectations, with one pre-existing
empty-test skip and 83 test diagnostic warnings. Test warning counts are
separate from check-level warnings. No diagnostic flags or package compiler
warnings were suppressed. Full test and compiler logs are retained.

## Exact submission archive

`genpca_0.2.2.tar.gz` (1,243,956 bytes) is the source archive built and checked
by the successful full R 4.6.1 check above, including the PDF manual and
rebuilt vignettes. It was downloaded from that run without rebuilding.

SHA-256:
`d06182849aa2f55dd51921924079823fbc309689a6f80defceb338fe3b45936d`

The artifact's checksum and 128 included source files were verified after
download. The local full check of this exact archive with current Matrix 1.7-6
also passed, with the three notes explained below.

Environment-specific builds have different bytes because of generated
vignette output and build metadata. The earlier supplemental local check
used its own archive, SHA-256:
`4f5e400e6a49724fc84024547cf13838adc8c1479cb3b249a121ff5361a60c64`.

## Local notes

1. Incoming feasibility identifies a new submission of an archived package.
   This resubmission addresses the failure that led to archival.
2. `unable to verify current time`: local network/time-service limitation.
3. HTML validation was skipped because the Mac's HTML Tidy is too old.
   The fresh hosted full R-release check passed HTML manual validation.

The first local check stopped during incoming feasibility when CRAN's
`PACKAGES.in` request timed out. The same local archive passed the complete
check on retry with a 180-second internet timeout; all diagnostic checks
remained enabled. Both logs were retained.

## Historical validation

PR #52 also passed regression and compiler/R-devel/R-release/manual checks
on 2026-10-05. The final historical workflow was
<https://github.com/bbuchsbaum/genpca/actions/runs/37345083392>.
Those results are historical context, separate from the fresh checks above.
The original reproduction and warning/environment audit are recorded in
<https://github.com/bbuchsbaum/genpca/pull/52>.

## Downstream dependencies

The CRAN package index retrieved on 2026-10-26 lists no reverse dependencies
on genpca in Depends, Imports, LinkingTo, Suggests or Enhances.
