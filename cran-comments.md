## Submission type

This is the first CRAN submission of `genpca`.

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Brad Buchsbaum <brad.buchsbaum@gmail.com>'
  New submission

  Expected for a first submission.

* checking C++ specification ... NOTE
  Specified C++14: please drop specification unless essential

  The C++14 specification is essential. The package's C++ sources use
  RcppArmadillo, which requires a C++14 compiler in current releases.
  `SystemRequirements: C++14` is declared in DESCRIPTION and
  `CXX_STD = CXX14` is set in `src/Makevars[.win]`.

## Test environments

* local macOS 14.2, R 4.5.3 (release)
* GitHub Actions ubuntu-latest, R 4.5 (release) -- pkgdown workflow
* win-builder R-devel and R-release (pending)
* macOS R-release builder (pending)

## Downstream dependencies

There are currently no reverse dependencies.

## Notes for the reviewer

* The package wraps published methods (Allen, Grosenick & Taylor
  2014; Abdi 2007). DOIs and references appear in the Description
  field as required.
* Linking dependencies (`RcppArmadillo`, `RcppEigen`, `RSpectra`)
  are declared under `LinkingTo`. The package explicitly requests
  C++14 for compatibility with current RcppArmadillo releases.
