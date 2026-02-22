# genpca: Vignette Guide

This is a short map of the genpca documentation. Each vignette is small
and focused; pick the one you need.

## Basics

[`vignette("gpca-basics", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-basics.md)
covers a quick start, identity vs diagonal metrics, choosing `ncomp`,
and a tiny real-data example (USArrests).

## Metrics

[`vignette("gpca-metrics", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-metrics.md)
shows how to build practical row/column metrics (AR(1), spatial kernels,
graph Laplacians, heteroscedastic weights), how
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
can learn metrics, and what to watch for with SPD remedies.

## Scale & Special Cases

[`vignette("gpca-scale", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gpca-scale.md)
gives backend choices (`eigen`, `spectra`, `deflation`), sparse
workflows, the covariance-only route (`genpca_cov`), out-of-sample
projection tips, and performance notes.

## PLS

For generalized PLS/CCA, see
[`vignette("gplssvd-reference", package = "genpca")`](https://bbuchsbaum.github.io/genpca/articles/gplssvd-reference.md)
for the mathematical reference and a quick-start code chunk.
