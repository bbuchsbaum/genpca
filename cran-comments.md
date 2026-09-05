# cran-comments

## Submission type

This is the first CRAN submission of `genpca` (version 0.2.0).

## R CMD check results

0 errors | 0 warnings | 2 notes

The same source tarball was checked with `R CMD check --as-cran`, including
vignettes and the PDF manual, against eigencore 1.0.3 and 1.3.0. Both runs
passed tests (1,403 expectations, zero failures), examples, vignette rebuilds
and PDF manual generation.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Brad Buchsbaum <brad.buchsbaum@gmail.com>'

  New submission

  Expected for a first submission.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
  Please obtain a recent version of HTML Tidy by downloading a binary
  release or compiling the source code from <https://www.html-tidy.org/>.

  This is a local tool limitation. The previous HTML math errors have been
  corrected; structural HTML validation with a current Tidy remains pending.

## Test environments

* Local macOS 14.3 (Apple Silicon), R 4.5.1, Apple clang 15.0.0,
  with CRAN eigencore 1.0.3 in an isolated library.
* Same environment with installed development eigencore 1.3.0.

A temporary Makevars selected Apple clang and the existing Homebrew Fortran
runtime. No compiler warnings were suppressed and user compiler settings were
unchanged. The ordinary Homebrew-clang configuration still emits the known
warning from R's own Boolean.h.

Not yet rerun for this candidate: Linux R-release/R-devel, Windows
R-release/R-devel, win-builder and macOS builder. Older results do not certify
this candidate; obtain and record these results before submission.

## Downstream dependencies

There are currently no CRAN reverse dependencies (new package).

## Notes for the reviewer

* The package implements published methods (Allen, Grosenick & Taylor, 2014,
  <doi:10.1080/01621459.2013.852978>; Abdi, 2007). References appear in the
  Description field.
* Compiled code uses RcppArmadillo and RcppEigen headers, declared under
  `LinkingTo`. Iterative R-level solves use eigencore (>= 1.0.3). RSpectra is
  no longer a dependency. The package uses R's default C++ standard and sets
  no `CXX_STD`.

## Candidate receipt (2026-09-04)

Source tarball SHA-256:
`3a73c06123b83fb7a61e4df81f77d46150b47ce54e049c4d39a558d33a910335`

Detailed logs: `plans/release-review-2026-09-04/fixed/` (build-ignored).
