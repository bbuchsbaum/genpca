# genpca audit remediation plan (2026-09-04)

Source: the "Section 8/9" findings from the Python-port design document. Each
item was checked against current `master` (1a0ff8d) by reading the code and,
where the claim is about runtime behaviour, by running it. A fresh-context
reviewer re-ran every script and re-checked every line reference; its
corrections are folded in below. Scripts live in the session scratchpad
(`t81.R` .. `t88.R`); their essential cases become regression tests in
Phase 0.

## 0. Execution status (2026-09-04)

Executed in the working tree (uncommitted): Phases 0-5 and 7 complete, the
eigencore swap folded into Phase 1's spectra work, Phase 6 done at the
documentation level with the pre-factored-metric input deferred. Two
fresh-context reviews (Phases 1-2 + swap; Phases 3-5) were applied. See
`PROGRESS.md` for the per-phase log and deviations.

## 1. Verdicts

| # | Claim | Verdict | Evidence |
|---|-------|---------|----------|
| 8.1 | `maxeig` truncates the metric before looking at X | **Confirmed** | `gmdLA()` → `compute_sqrtm()` uses `RSpectra::eigs_sym(M, k = maxeig, "LM")` and builds `sqrtm`/`invsqrtm` from the retained eigenpairs only (`R/gpca.R:946-985`); only the small-side metric is factored (R in primal, Q in dual). `maxeig` is used nowhere else in R/ or src/. Runtime (t81): p = 6, A with eigenvalues (1, .9, .8, .7, .6, 1e-2), X with variance 1e4 along the small-eigenvalue direction. `maxeig = 800` gives sdev 63.56, 6.97, 6.64 (matches the exact eigenvalues); `maxeig = 4` gives 6.86, 6.54, 6.27 and a PC1 with \|cos\| = 0.04 to the true PC1. The dominant component is gone; only a "result is approximate" warning is emitted. |
| 8.2 | PSD/SPD conflated; `tol` not forwarded; vector path admits `-tol` | **Confirmed** | `is_spd()` is documented and implemented as a PSD-within-tolerance probe (Cholesky of `A + tol*scale*I`, `constraints_utils.R:15-44`). `ensure_spd()` calls `is_spd(M)` three times with the default `tol`, never its own argument (`constraints_utils.R:118,136,148`). Runtime (t82): `is_spd(0_3x3)` is TRUE and `ensure_spd(0_3x3)` returns the zero matrix; `ensure_spd(diag(1,-5e-7), tol = 1e-12)` returns the input with eigenvalue -5e-7. Vector path: `prep_constraints()` accepts `A >= -tol` (1e-6, `gpca.R:14,51`) and puts the value on the diagonal; `genpca()` never passes a tolerance to `gmdLA()` (`gpca.R:499-505`), whose diagonal cutoff is 1e-8 (`:927`). With `A = c(1,1,1,1,-5e-7)`, `method = "eigen"` errors ("R (diagonal) must be PSD") while `spectra`, `deflation`, `randomized` run silently. Backend-dependent behaviour, as claimed. |
| 8.3 | Silent one-triangle symmetrization; tolerant `clip_psd` | **Confirmed** | `forceSymmetric(., uplo = "U")` at `constraints_utils.R:82,114` and `internal_ops.R:83`; `genpca_cov*()` use default `forceSymmetric()` (`gpca_cov.R:171,182,275,302`). Only `remedy = "error"`/`"identity"` check symmetry first (`gpca.R:23-37`); the default `"ridge"` does not. Runtime (t83): `A = [[2, .1],[.9, 2]]` → `genpca()` gives sdev 9.21 = upper-triangle answer (lower: 10.81, averaged: 10.04) with zero warnings/messages. Same for `genpca_cov()` (6.436 vs 6.359). `clip_psd(diag(1,-5e-7))` returns min eigenvalue -5e-7 because its fast path is `is_spd(M, tol)`. |
| 8.4 | Absolute tolerances break scale equivariance | **Confirmed, and worse than stated** | `gmdLA(tol = 1e-8)` applies that absolute value to metric eigenvalues (`:927`), target eigenvalues (`:1025,1073`), `norms_sq > tol^2` (`:1045,1094`) and `total_variance < tol` (`:1113`). `genpca()` uses `total_variance < 1e-8` for spectra/randomized (`:521,577`) and `> 1e-8` in the deflation fallback (`:484`). Inside "spectra" there are two different cutoffs: the diagonal-metric R path keeps `d > tol_spectra` (1e-9, `gmd_fast.R:55`), the general-metric C++ path filters eigenvalues, `d² > tol` i.e. `d > 3.2e-5` (`gmd_fast.cpp:223,237`), and the same `tol` is also the Spectra convergence tolerance (`:106`). Runtime (t84, 60x8, general A, diagonal M, ncomp 4, c = 1, 1e-3, 1e-5, 1e-7): eigen finds 4/4/1/error; spectra 4/4/4/0; deflation 4/4/4/4; randomized 4/≈/wrong/0. Randomized is already ~2% off at c = 1e-3 (sdev/c 15.30, 9.78, 8.35, 6.84 vs 15.32, 9.83, 8.44, 6.98, no warning) and at c = 1e-5 returns 0.036, 0.0095, 0.0060, 0.0034: the absolute `jitter_metric = 1e-10` added to the Gram matrix in `metric_orthonormalize()` swamps a Gram that scales like c^(2(2q+1)) under power iteration. Setting `jitter_metric = 1e-20` restores the correct answer (t88). Deflation is already relative (NEWS 0.1.0; `gpca.cpp:83,137`) and is the only scale-equivariant backend. |
| 8.5 | Covariance interface mixes GMD and generalized eigen | **Confirmed as description; redesign judgment** | `genpca_cov(method = c("gmd","geigen"))` dispatches to two different problems (eigen of `R^{1/2} C R^{1/2}` vs `C v = λ R v`). The GMD branch force-symmetrizes `C` and `R`, clamps metric eigenvalues with `pmax(., 0)` regardless of magnitude (diagonal path `gpca_cov.R:190`, general `:198`), clamps target eigenvalues the same way (`:223`) and never checks that `C` is PSD. The docs do say the two methods differ. Whether to split them is an API decision, but the validation gaps are real. |
| 8.6 | Sparse metrics densified before Cholesky | **Confirmed for the spectra path; not for deflation** | `gmd_fast_cpp()` routes a sparse small-side metric to `gmd_fast_cpp_sp()` (`gmd_fast.R:184-196,201-212`), whose `gmd_fast_auto()` does `arma::mat Rdense(R); arma::chol(L, Rdense)` (`gmd_fast.cpp:490-497`); the cached path `get_chol_lower_dense()` also densifies (`gmd_cache.R:49`). Moreover, the spectra and randomized branches densify X itself before calling C++ (`gpca.R:511,555`), so under those methods nothing is sparse end-to-end, diagonal metric or not. The C++ deflation kernel (`gpca.cpp:34-160`) takes `sp_mat` X, Q, R and only multiplies, no factorization anywhere in the file, so deflation is the one sparse end-to-end path. Validation itself runs one CHOLMOD factorization of a sparse metric under every remedy (`constraints_utils.R:34`), which costs fill-in, not densification. |
| 8.7 | `gpca_mle` exit rescale is inconsistent with its penalized objective | **Confirmed** | Algebra: `M ← sM, A ← A/s` leaves `p·logdet Σ_r + n·logdet Σ_c + tr(M E A E')` unchanged but turns `pλ tr M + nλ tr A` into `pλs tr M + nλ tr A / s`. The rescale is at `gpca.R:870-882`; `loglik = tail(loglik_path, 1)` at `:896` is never recomputed. Runtime (t87, 30x6, ncomp 2, λ = 1e-3, 8 iterations): `scale_fix = "trace"` and `"none"` report the identical loglik 880.55, but "trace" rescales by s = 509; the penalty goes from 312 to 79 420, i.e. the objective at the returned (M, A) is about 39 554 lower than reported. The penalty-optimal scale given the shapes is s* = 0.9994, so the flip-flop iterate was already at the right scale and the default `"trace"` moves it away. |
| 8.8 | Default `"ridge"` silently alters the supplied metric | **Confirmed** | `genpca(constraints_remedy = c("ridge", ...))`. The ridge branch (`gpca.R:41-42,78-79`) and `ensure_spd()` emit nothing at any verbosity; `genpca(verbose = TRUE)` prints only "Preparing constraints...". Runtime (t88): `A = [[1,2],[2,1]]` (eigenvalues -1, 3) is replaced by `[[2.1,2],[2,2.1]]` with zero warnings or messages. The same silent repair reaches the PLS family with no opt-out at all: `.metric_operators()` calls `ensure_spd(W)` unconditionally for any general metric (`internal_ops.R:84`) and silently clamps negative diagonal weights (`:62`); it is used by `genpls()`, `genplsc()`, `gplssvd_op()` and `rpls()` via `plsutils.R:121,146`, `gplssvd_op.R:131-137`, `weight_operators.R:40`. |

Section 9 historical items: NEWS 0.1.0 records fixes for scores = U D,
`reconstruct(colind =)`, relative deflation thresholds, real spectral clip,
dense-PSD validation, `.Random.seed` hygiene, and exact-byte cache keys. Code
matches (`gpca.R:624-626`, `gpca.R:1439-1443`, `gmd_cache.R:11,28-40`).
These stay as regression tests; none is a current defect.

## 2. Remediation phases

Ordering: tests first, then the numerical-semantics core (Phases 1-2, which
change results), then API/policy changes (Phases 3-5), then the performance
item (Phase 6). Each phase ends with `devtools::test()` and, for phases that
change exported behaviour, `R CMD check`. Version bumps to 0.2.0 at Phase 3.

### Phase 0. Regression tests that fail today

File: `tests/testthat/test-audit-2026-09.R`. One `test_that()` per row above,
written so it fails on current master and passes once the corresponding phase
lands (`skip()` per phase if the suite must stay green in between).

- 8.1: the 6-column counterexample; `genpca(maxeig = 4)` must match
  `maxeig = 800` to 1e-6 in sdev (after Phase 2 that means: error unless
  `maxeig = Inf`, and `maxeig = Inf` matches).
- 8.2: `ensure_spd(zero)` must return a PD matrix; vector weights of -5e-7
  must be treated identically by all four backends.
- 8.3: asymmetric `A` with relative asymmetry 0.2 must error under every
  remedy; asymmetry 1e-13 must be averaged; `clip_psd(diag(1,-5e-7))` must have
  min eigenvalue ≥ 0.
- 8.4: for each backend separately, with a fixed seed, `genpca(c*X)` for c in
  {1, 1e-3, 1e-5, 1e-7} must return the same number of components and sdev/c
  within 1e-6 of the c = 1 result. Cross-backend agreement is a different,
  looser test (`test-randomized-backend.R` uses a variance ratio > 0.9).
- 8.7: `gpca_mle()$loglik` must equal the penalized objective evaluated at the
  returned `(M, A, fit)` for every `scale_fix`.
- 8.8: remediation must be observable: a warning of a catchable class
  carrying the diagnostic report.

### Phase 1. PSD/PD semantics, symmetry, and relative tolerances

Files: `R/constraints_utils.R`, `R/gpca.R` (`prep_constraints`, `gmdLA`,
`genpca` variance checks), `R/internal_ops.R`, `R/gmd_fast.R`
(`gmd_fast_cpp`, `gmd_fast_diag`, `metric_orthonormalize`),
`src/gmd_fast.cpp` (separate `tol` from `rank_rtol`), `R/gpca_cov.R`.

1. Split `is_spd()` into `is_psd(A, rtol)` and `is_pd(A, rtol)`. Both stay
   shifted-Cholesky probes so large sparse input never needs an
   eigendecomposition: `is_psd` factors `A + rtol·scale·I` (λ_min > -rtol·scale),
   `is_pd` factors `A - rtol·scale·I` (λ_min > +rtol·scale). Keep `is_spd` as
   a deprecated alias of `is_psd` for one release.
2. `ensure_spd()` forwards `tol` everywhere it probes, and uses `is_pd` for its
   "already fine" exit. Any λ_min ≤ rtol·scale therefore gets the shift, so a
   near-singular PSD input is shifted deterministically rather than flipping
   with roundoff (a matrix within roundoff of the boundary receives a shift
   of order rtol·scale either way).
3. `.metric_operators()` stops calling `ensure_spd()`: it symmetrizes per
   step 4, checks `is_psd` and errors on failure, and keeps its eigen
   pseudo-inverse for PSD-singular weights (so Laplacian-type PLS weights keep
   today's behaviour instead of acquiring a ridge shift from step 2).
   Remedy plumbing for the PLS family comes in Phase 3.
4. Vector/diagonal weights: accept exactly nonnegative entries; entries in
   `[-rtol·max|w|, 0)` are set to exactly 0 with a message; anything below is
   an error. `rtol` default `sqrt(.Machine$double.eps)`. Remove the divergent
   backend-specific handling: the `gmdLA` 1e-8 check (`gpca.R:927`), and the
   silent `pmax(., 0)` clamps in `gmd_fast_diag()` (`gmd_fast.R:26-27`),
   `.metric_operators()` (`internal_ops.R:62`) and `genpca_cov_gmd()`
   (`gpca_cov.R:190`).
5. Symmetry: new helper `symmetrize_or_stop(A, rtol, name)` computing
   `‖A - Aᵀ‖_F / max(‖A‖_F, eps)`; average the triangles below `rtol`
   (default 1e-10), else stop with the measured asymmetry. Keep the existing
   message prefix "Matrix A must be symmetric" (matched by
   `test-coverage-gpca-branches.R:6,10`). Use it in `prep_constraints()`
   (all remedies, before any repair), `ensure_spd()`, `clip_psd()`,
   `.metric_operators()`, both `genpca_cov*()` branches. Asymmetry is an
   input error, not something a PSD remedy should paper over.
   Test impact: `test-coverage-gpca-branches.R:214-222` feeds `genpca_cov()`
   a matrix with relative asymmetry 0.046 and expects success; rewrite it to
   expect the error and add a 1e-13-asymmetry case that succeeds.
6. `clip_psd()`: drop the tolerant fast path; always eigendecompose when
   below `dense_maxn`. Contract: output has min eigenvalue ≥ 0 exactly.
7. Rank decisions relative to the leading value:
   - `gmdLA`: metric eigenvalues kept if `> rtol · max(values)`; target
     eigenvalues `> rtol · λ_1`; `norms_sq` floor relative to
     `total_variance`; total-variance test becomes `!is.finite || <= 0`.
   - `genpca()`: replace all three absolute total-variance tests
     (`:484,521,577`) with the same finite/positive test; component filter
     `d > rtol · d_1` applied uniformly after every backend.
   - C++ `gmd_fast.cpp`: add a `rank_rtol` argument; keep `tol` only as the
     Spectra convergence tolerance; `sqrt_pos()` and the `ou`/`ov` division
     guards use `rank_rtol · λ_1`. Expose `rank_rtol` from `gmd_fast_cpp()`
     and use it in `gmd_fast_diag()` too, so "spectra" has one cutoff.
   - `metric_orthonormalize()` (R and C++): `jitter` becomes relative,
     `jitter · max(diag(G))`, so `jitter_metric = 1e-10` means 1e-10 of the
     Gram scale. This fixes the randomized wrong-answer case.
   - `genpca_cov*()`: `tol` → `rank_rtol` semantics for `R` eigenvalues and
     target eigenvalues; validation replaces clamping (error below
     `-rtol·λ_1`, set to 0 above).
8. Document the convention once (`?genpca` "Tolerances" section): absolute
   tolerances govern iteration, relative tolerances govern rank.

Verification: Phase 0 tests for 8.2/8.3/8.4; `test-method-differential.R`,
`test-deflation-vs-eigen.R`, `test_prep_constraints.R`, `test-ensure-spd.R`,
`test-coverage-gpca-branches.R` must pass, updating only assertions that
encoded the old absolute cutoffs or the one-triangle behaviour.

### Phase 2. Remove metric-first truncation (8.1)

File: `R/gpca.R` (`gmdLA`, `compute_sqrtm`, the `method = "auto"` heuristic),
docs for `maxeig`/`warn_approx`.

1. Delete the `eigs_sym(k = maxeig)` branch. For a general metric the exact
   options are: (a) `is_pd` → Cholesky `R = L Lᵀ`, target `Lᵀ (XᵀQX) L`
   (similar to `(XᵀQX) R`, hence the same eigenvalues as
   `R^{1/2} XᵀQX R^{1/2}`), map back `V = L^{-ᵀ} Z` so `VᵀRV = ZᵀZ = I`;
   (b) PSD-singular → full symmetric eigendecomposition of the metric (O(p³)
   but exact). This is what `gmd_primal_impl` already does
   (`gmd_fast.cpp:225`); reuse `get_chol_lower_dense()` for (a).
2. `maxeig` becomes a cost guard: if the small-side metric is general and its
   dimension exceeds `maxeig`, `method = "eigen"` stops with a message
   recommending `method = "spectra"` or `"deflation"`, unless the user passes
   `maxeig = Inf`. Never approximate silently. Deprecate `warn_approx` with a
   one-release warning.
3. `method = "auto"` must route past the guard: today it falls back to
   `"eigen"` for general metrics with `min_dim < 1000` (`gpca.R:405-411`), so
   800 < min_dim < 1000 with a dense general metric would newly error. Add
   "small-side metric is general and its dimension > maxeig ⇒ spectra" to the
   heuristic.
4. `test-method-differential.R:80-108` exists to exercise the RSpectra
   truncation branch (dual, Q = 140 > `maxeig = 80`); rewrite it to assert
   (i) `eigen` errors at `maxeig = 80`, (ii) `eigen` with `maxeig = Inf`
   matches `spectra` to 1e-6.
5. Cache the factor on the metric attribute as today.

Verification: Phase 0 test for 8.1; `test-coverage-gpca-branches.R` and any
test passing `maxeig`/`warn_approx` updated; eigen-vs-spectra equivalence at
p > maxeig with a dense SPD metric.

### Phase 3. Explicit metric repair (8.8)

Files: `R/gpca.R` (`genpca`, `prep_constraints`, `gpca_mle`),
`R/constraints_utils.R`, new `R/repair_metric.R`, `R/internal_ops.R`,
`R/plsutils.R`, `R/weight_operators.R`, `R/genpls.R`, `R/genplsc.R`,
`R/gplssvd_op.R`, `R/rpls.R`, NEWS, vignette `gpca-metrics`.

1. Default `constraints_remedy = "error"` (breaking; NEWS "behaviour change";
   bump to 0.2.0). PSD-singular metrics remain valid input under "error"
   (GMD is defined for PSD); only indefinite or asymmetric input is rejected.
   Test impact: `test_gpca.R:148-153` uses `adjoin::temporal_adjacency()`
   as M, which is indefinite (min eigenvalue -0.21 at n = 200) and relies on
   the silent ridge; pass `constraints_remedy = "ridge"` explicitly there and
   audit every other test/example/vignette that supplies an adjacency matrix.
2. Export `repair_metric(A, method = c("ridge", "clip"), rtol)` returning the
   repaired `Matrix` with attribute `"repair_report"`: original min
   eigenvalue (or Gershgorin bound when too large to eigendecompose), applied
   shift, resulting condition number (when computable), rank. Print method.
3. When a non-"error" remedy *changes* the matrix inside `genpca()`, emit a
   warning of class `genpca_metric_repaired` carrying the report; no output
   when the input was already valid.
4. `gpca_mle()` drops its hard-coded `constraints_remedy = "ridge"`
   (`gpca.R:790`) and inherits the new default; its metrics are inverses of
   SPD matrices, so no repair fires and no warning appears.
5. PLS family: `genpls()`, `genplsc()`, `gplssvd_op()`, `rpls()` gain the
   same `constraints_remedy` argument (default `"error"`), threaded through
   `plsutils.R` and `weight_operators.R` into `.metric_operators()`, which
   applies `repair_metric()` only when asked.
6. Update every example/vignette that relies on the silent ridge.

Verification: Phase 0 test for 8.8; full suite; `R CMD check` examples.

### Phase 4. Split the covariance interface (8.5)

Files: `R/gpca_cov.R`, NAMESPACE, `_pkgdown.yml`, tests
`test_genpca_cov*.R`, `test-coverage-gpca-branches.R:229-270`.

1. `genpca_cov(C, R, ncomp, ...)` becomes GMD-only. It requires `C` PSD
   within `rtol` (it models `C = XᵀMX`); indefinite `C` is an error.
2. New export `geigen_cov(C, R, ncomp, ...)` for `C v = λ R v`, documented
   as a different estimator with its own optimality statement. The
   generalized eigenproblem is well defined for indefinite `C`, so it keeps
   today's warning-and-proceed behaviour for non-PSD `C`
   (`test-coverage-gpca-branches.R:229-237` expects that warning), and keeps
   its `constraints_remedy` for `R` (`test_genpca_cov.R:127-139`,
   `test-coverage-gpca-branches.R:240-270`).
3. `genpca_cov(method = "geigen")` kept for one release, emitting a
   deprecation warning and forwarding to `geigen_cov()`.
4. Both apply the Phase 1 symmetry threshold and relative rank tolerance.

Verification: `test_genpca_cov_equivalence.R`, `test_genpca_cov_gmd.R`,
`test_genpca_cov.R`, `test-coverage-gpca-branches.R` rerun against the split
functions, with the geigen remedy tests re-pointed at `geigen_cov()`.

### Phase 5. `gpca_mle` objective consistency (8.7)

File: `R/gpca.R` (`gpca_mle`), `tests/testthat/test-gpca_mle.R`,
`tests/testthat/test-phase3-hardening.R`.

Decision: keep the function (it is exported and documented) but make it honest.

1. Default `scale_fix = "none"`. The penalty already identifies the scale
   (empirically s* ≈ 1 at convergence), so the flip-flop output is the
   penalized optimum; "trace"/"det" are post-hoc reparameterizations.
2. `loglik` is always the penalized objective at the returned `(M, A, fit)`.
   `loglik_path` stays the per-iteration path (monotone). When
   `scale_fix != "none"`, also return `loglik_unpenalized` and
   `loglik_rescale_delta`, and document that the penalized objective is
   *not* invariant to the rescale.
   Test impact: `test-phase3-hardening.R:20` asserts
   `loglik == tail(loglik_path, 1)` for all three `scale_fix` values; keep it
   for `"none"`, and for `"trace"`/`"det"` assert `loglik` equals the
   objective recomputed at the returned metrics instead.
3. Add `experimental` to the Rd title and a "not identifiable with
   `lambda = 0`" note; keep the `lambda = 0` warning.

Verification: Phase 0 test for 8.7; monotone `loglik_path` test retained.

### Phase 6. Sparse metrics without densification (8.6)

Files: `R/gmd_fast.R`, `src/gmd_fast.cpp`, `R/gmd_cache.R`, `R/gpca.R`
(spectra/randomized branches), docs (`gpca-scale` vignette).

1. Documentation first: state precisely what is sparse end-to-end today
   (only `method = "deflation"`, which keeps X, Q, R sparse in C++), that
   `"spectra"` and `"randomized"` densify X before calling C++
   (`gpca.R:511,555`) and factor a general metric densely, and that
   validation performs one sparse CHOLMOD factorization of a sparse metric
   under every remedy (fill-in, not densification).
2. Add an explicit `metric_factor` input path (class `"metric_factor"`
   wrapping a lower-triangular `L` with `A = L Lᵀ`, dense or `dtCMatrix`), so
   users who already hold a sparse Cholesky can pass it.
3. Sparse CHOLMOD path: when the small-side metric is sparse and `is_pd`,
   factor it with `Matrix::Cholesky(perm = TRUE)` and pass the permuted
   sparse `L` to a new `gmd_primal_impl` overload templated on `arma::sp_mat`
   for `L`; the operator applies `L`/`Lᵀ` via sparse triangular solves.
   Densify only when the factor's fill-in exceeds a threshold, with a
   message. Let the spectra branch pass a sparse X through
   (`gmd_fast_cpp_*` overloads on `sp_mat` X exist for deflation; add them
   for spectra) so the promise "sparse X + diagonal or factored metric ⇒
   sparse" holds.
4. `get_chol_lower_dense()` stops silently densifying sparse input; callers
   choose.

Verification: memory test with `bench::mark(memory = TRUE)` on a
20 000-column banded metric; equivalence against the dense path on a small
problem.

### Phase 7. Documentation and naming

- `sdev` normalization: documentation-only per project memory (the convention
  is set at the `multivarious` level). Add one sentence to `?genpca` "Value"
  stating `sdev` = GMD singular values, not standard deviations, with the
  divisor a user would apply.
- NEWS entries for every phase; version 0.2.0 at Phase 3.
- Refresh `gpca-metrics` and `gpca-scale` vignettes for the new defaults.

## 3. Deliberately not adopted from the Python decisions

- Removing `gpca_mle`: it is exported and documented; Phase 5 makes it
  consistent instead.
- Dropping the `method` switch entirely: `"auto"` heuristic stays; only the
  covariance API is split.
- Global cache removal: the Cholesky cache is bounded (LRU 16) and keyed on
  exact bytes; leave it, but Phase 6 stops it from densifying.
- "Metric-aware centering", "primal/dual feature components exposed
  separately", single common row metric in PLS: API redesigns for the port,
  not defects in R; out of scope here.

## 4. Effort estimate

| Phase | Size | Risk |
|-------|------|------|
| 0 | small | none (tests only) |
| 1 | large | changes numerical results near rank boundaries; most tests touched |
| 2 | medium | exact by construction; p > maxeig with a general metric becomes a hard stop for `eigen` (was a warning) and `auto` must route around it |
| 3 | medium | breaking default; user-visible across genpca and the PLS family |
| 4 | small-medium | API churn, deprecation path |
| 5 | small | none beyond default change |
| 6 | large | C++ + CHOLMOD; can be deferred after step 1 (docs) |
| 7 | small | none |
