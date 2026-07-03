# genpca 0.1.0

* Initial CRAN submission.
* Generalized PCA with row metric `M` and column metric `A`
  (`genpca()`), following Allen, Grosenick & Taylor (2014).
* Covariance-based GPCA from a precomputed `C = X' M X`
  (`genpca_cov()`), with `eigen` and `gmd` paths.
* Multiple computational backends: `eigen`, `spectra` (matrix-free C++
  via RSpectra), `randomized`, and `deflation`. The `auto` heuristic
  picks among them.
* Generalized PLS / PLS-SVD on two blocks (`genpls()`, `genplsc()`)
  and an operator-level interface (`gplssvd_op()`) that avoids
  materialising whitened matrices.
* Sparse functional PCA (`sfpca()`), regularised PLS (`rpls()`), and
  matrix-normal PCA via maximum residual likelihood (`mnpca_mrl()`).
  `sfpca()` solves each rank-1 subproblem exactly via C++ coordinate
  descent in the constraint form of Allen & Weylandt (2019), deflates
  implicitly (sparse inputs never densify), and selects sparsity
  penalties per component by BIC along a warm-started `lambda` path;
  smoothness weights default to the scale-free `1 / lambda_max(Omega)`.
  The former `uthresh`/`vthresh` quantile heuristics are deprecated.
  `sfpca()` now returns a `multivarious` `bi_projector` (class `"sfpca"`),
  so `scores()`, `components()`, `sdev()`, and `reconstruct()` work as
  they do for `genpca()`. The pre-0.1 list fields `$d` and `$u` remain
  readable but emit a deprecation warning.
* Maximum-likelihood metric learning (`gpca_mle()`).
* SPD constraint handling with `"ridge"`, `"clip"`, and `"identity"`
  remediation strategies.
* Vignettes covering getting started, metric recipes, scaling, and a
  reference implementation for GPLSSVD.
