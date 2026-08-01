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

### Bug fixes and behavior changes

* `gpca_mle()` now optimizes a single penalized (MAP) matrix-normal
  objective with exact block updates (sequential flip-flop covariance
  updates; the ridge `lambda` is part of the objective), so
  `loglik_path` is monotone non-decreasing. `scale_fix` is applied once
  at exit as a joint, likelihood-preserving rescale (row covariance
  normalized, factor absorbed into the column covariance) instead of
  normalizing both factors independently during iteration, and the
  returned `fit` is computed with the returned canonicalized metrics.
  Reported log-likelihood values now include the penalty term, so they
  differ numerically from previous versions.
* `scores()` on `genpca` fits now matches Allen, Grosenick & Taylor (2014)'s
  convention (`z_k = X A ov_k = ou_k d_k`) and is identical to
  `project(fit, X)` on the training data; previously the two could disagree.
* `reconstruct.genpca(fit, colind = )` now applies the inverse
  pre-processing transform to the selected columns, rather than to the
  first `length(colind)` columns.
* Deflation convergence and degeneracy thresholds (the `threshold` argument
  of `genpca(method = "deflation")` and the internal SFPCA solver) are now
  relative to the scale of the problem rather than absolute, so results are
  invariant to rescaling the input data.
* `constraints_remedy = "clip"` now performs a real spectral clip to the
  PSD cone (previously it did not clip); it requires a dense
  eigendecomposition and refuses sparse input larger than 2000 rows/columns,
  recommending `"ridge"` instead.
* Fixed PSD validation for dense (non-diagonal) constraint matrices so that
  indefinite input is reliably rejected/repaired rather than silently
  accepted.
* `genpca(method = "randomized")` no longer mutates the caller's
  `.Random.seed`; the C++ kernel seeds its own stream from
  `seed_randomized`, which now fully controls reproducibility.
* Added input validation (dimensions, finiteness, penalty/`scad_a` range,
  strictly positive diagonal) to the internal C++ coordinate-descent solver
  used by `sfpca()`.
* Fixed a cache-key collision in the internal matrix-decomposition cache:
  the digest previously used rounded matrix entries, which could collide
  for distinct matrices (e.g. a metric and its ridge-remediated version)
  and return the wrong cached factor; the digest now hashes exact bytes.
* `genpls()`'s `$ncomp` field now reports the number of components actually
  extracted, which may be less than the requested `ncomp`.
