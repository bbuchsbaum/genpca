# genpca 0.2.2 check notes (candidate)

## Reported failure

CRAN's Fedora R-devel checks of 0.2.1 report `starting vector near the null
space` in `test-gplssvd-op-large.R`. irlba 2.4.1 removed the `mult` argument;
the old implementation supplied a zero placeholder matrix and a callback
through that argument. The replacement supplies a matrix-free S4 operator
with forward and adjoint multiplication methods using irlba's documented
interface. No tests were skipped, seeds changed, or tolerances loosened.

## Validation status

Candidate validation is in progress. This file does not claim successful
R-hub or full CRAN checks until their logs have been reviewed. The previous
0.2.1 submission results are historical and do not validate this candidate.
