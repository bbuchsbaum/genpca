# cran-comments

## Submission type

This is a resubmission of `genpca` (version 0.2.1). It addresses the test
failure reported on the CRAN Fedora/OpenBLAS check.

## Response to the previous check

The reported failure was:

```
Expected `genpca:::.mnpca_safe_solve_spd(matrix(NaN, 2, 2))` to throw an error.
```

The internal helper had relied on `chol()`/`chol2inv()` to signal an error for
non-finite input. That behavior differed across numerical backends. The helper
now rejects non-finite matrices before factorization and also refuses a
non-finite inverse. No test was skipped and no tolerance was loosened.

## R CMD check results

Hosted R-hub Fedora R-devel check: 0 errors | 0 warnings | 0 notes

Final `genpca_0.2.1.tar.gz`, local `R CMD check --as-cran`:
0 errors | 1 warning | 3 notes

The corrected code passed all tests on R-hub's `gcc16` platform: 1,588
expectations, zero failures and one intentionally empty test skipped. The
standard R-hub check used `--as-cran --no-manual --no-build-vignettes`.
The hosted run tested the exact 0.2.1 package-content commit `66053cf`:
<https://github.com/bbuchsbaum/genpca/actions/runs/35101331539>.
The final tarball passed all 1,588 expectations locally and rebuilt all
vignettes and the PDF manual.

* checking whether package `genpca` can be installed ... WARNING
  The local Homebrew clang reports `-Wfixed-enum-extension` as an unknown
  warning option while compiling R's own `R_ext/Boolean.h`. The warning does
  not originate in package code and was not suppressed.

* checking CRAN incoming feasibility ... NOTE
  `Days since last update: 1` is expected for this corrective release.

* checking for future file timestamps ... NOTE
  The local check was unable to verify current time. This is an environmental
  network/time-service limitation.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
  Please obtain a recent version of HTML Tidy by downloading a binary
  release or compiling the source code from <https://www.html-tidy.org/>.

  This is a local tool limitation. The previous HTML math errors have been
  corrected; structural HTML validation with a current Tidy remains pending.

## Test environments

* R-hub `gcc16`: Fedora Linux 44 (x86_64), R-devel r90540,
  GCC/GFortran 16.2.1, OpenBLAS 0.3.29: status OK.
* Local macOS 14.3 (Apple Silicon), R 4.5.1, Homebrew clang 20.1.8.

The local build used Homebrew clang and the existing Homebrew Fortran runtime.
No compiler warnings were suppressed and user compiler settings were
unchanged.

An earlier Windows R-devel check passed, as reported by the maintainer. This
followed fixes for randomized sketch rank loss and separation of the MLE
rescaling diagnostic from final refit/reevaluation differences. Both fixes
are covered by additional regressions. The successful Windows log is not
stored in the repository.

Additional environments not checked for this candidate: Linux R-release,
Windows R-release and the CRAN macOS builder.

## Downstream dependencies

There are currently no known CRAN reverse dependencies.

## Notes for the reviewer

* The package implements published methods (Allen, Grosenick & Taylor, 2014,
  <doi:10.1080/01621459.2013.852978>; Abdi, 2007). References appear in the
  Description field. `Abdi` and `Grosenick`, flagged by the Windows incoming
  spelling check, are the surnames of the cited authors.
* Compiled code uses RcppArmadillo and RcppEigen headers, declared under
  `LinkingTo`. Iterative R-level solves use eigencore (>= 1.0.3). RSpectra is
  no longer a dependency. The package uses R's default C++ standard and sets
  no `CXX_STD`.

## Candidate receipt (2026-09-16)

Source tarball SHA-256:
`da73e5b7aa7a919bb84babc7cc9aaa25f2c9cf6bccd0bd1289b24d679614511f`
