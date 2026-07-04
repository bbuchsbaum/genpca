# Changelog

## genpca 0.1.0

- Initial CRAN submission.
- Generalized PCA with row metric `M` and column metric `A`
  ([`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)),
  following Allen, Grosenick & Taylor (2014).
- Covariance-based GPCA from a precomputed `C = X' M X`
  ([`genpca_cov()`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)),
  with `eigen` and `gmd` paths.
- Multiple computational backends: `eigen`, `spectra` (matrix-free C++
  via RSpectra), `randomized`, and `deflation`. The `auto` heuristic
  picks among them.
- Generalized PLS / PLS-SVD on two blocks
  ([`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md),
  [`genplsc()`](https://bbuchsbaum.github.io/genpca/reference/genplsc.md))
  and an operator-level interface
  ([`gplssvd_op()`](https://bbuchsbaum.github.io/genpca/reference/gplssvd_op.md))
  that avoids materialising whitened matrices.
- Sparse functional PCA
  ([`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)),
  regularised PLS
  ([`rpls()`](https://bbuchsbaum.github.io/genpca/reference/rpls.md)),
  and matrix-normal PCA via maximum residual likelihood
  ([`mnpca_mrl()`](https://bbuchsbaum.github.io/genpca/reference/mnpca_mrl.md)).
  [`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
  solves each rank-1 subproblem exactly via C++ coordinate descent in
  the constraint form of Allen & Weylandt (2019), deflates implicitly
  (sparse inputs never densify), and selects sparsity penalties per
  component by BIC along a warm-started `lambda` path; smoothness
  weights default to the scale-free `1 / lambda_max(Omega)`. The former
  `uthresh`/`vthresh` quantile heuristics are deprecated.
  [`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
  now returns a `multivarious` `bi_projector` (class `"sfpca"`), so
  `scores()`,
  [`components()`](https://bbuchsbaum.github.io/multivarious/reference/components.html),
  `sdev()`, and
  [`reconstruct()`](https://bbuchsbaum.github.io/multivarious/reference/reconstruct.html)
  work as they do for
  [`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md).
  The pre-0.1 list fields `$d` and `$u` remain readable but emit a
  deprecation warning.
- Maximum-likelihood metric learning
  ([`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)).
- SPD constraint handling with `"ridge"`, `"clip"`, and `"identity"`
  remediation strategies.
- Vignettes covering getting started, metric recipes, scaling, and a
  reference implementation for GPLSSVD.
