# cran-comments

## Submission type

This is the first CRAN submission of `genpca`.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Brad Buchsbaum <brad.buchsbaum@gmail.com>'
  New submission

  Expected for a first submission.

<!-- TODO(Phase 5): re-run rcmdcheck --as-cran on the release tarball and
     update the counts/notes above verbatim before submitting. Local runs on
     a Homebrew-clang toolchain show a spurious install-time warning
     (-Wunknown-warning-option from R's own headers) that does not occur
     under Apple clang or gcc; verify it is absent on win-builder/mac-builder
     and do not mention it here if so. -->

## Test environments

* local macOS 14.3 (Apple Silicon), R 4.5.1
* win-builder R-devel and R-release
* macOS R-release builder (mac.r-project.org)

## Downstream dependencies

There are currently no reverse dependencies.

## Notes for the reviewer

* The package implements published methods (Allen, Grosenick & Taylor, 2014,
  <doi:10.1080/01621459.2013.852978>; Abdi, 2007). References appear in the
  Description field.
* Compiled code links against RcppArmadillo, RcppEigen, and RSpectra headers
  (declared under `LinkingTo`); the package builds with the default C++
  standard and sets no `CXX_STD`.
