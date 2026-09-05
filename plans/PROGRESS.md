# Execution log: audit remediation + eigencore swap

Plans: `2026-09-04-audit-remediation.md`, `2026-09-04-eigencore-swap.md`.
Sequencing decided 2026-09-04: Phase 0 → Phase 1 (with the eigencore swap
folded into the spectra-path work so `src/gmd_fast.cpp`'s Spectra kernel is
touched once) → Phase 2 → review → Phases 3, 4, 5 → Phase 6 → Phase 7 → R CMD check.
No commits unless asked; working tree carries the changes.

| Phase | Status | Notes |
|---|---|---|
| 0 tests | done | `tests/testthat/test-audit-2026-09.R`; 8.1-8.4 active, 8.7/8.8 skipped until phases 5/3 |
| 1 tolerances/PSD | done | is_psd/is_pd split, symmetrize_or_stop, strict clip, .clamp_weights, relative rank_rtol/metric_rtol everywhere incl. C++; suite green (2 pre-existing rpls failures only) |
| swap (1a/1b/1c) | done | R/solver_backend.R wraps eigencore; spectra path = SVD of whitened operator (R/gmd_fast.R .metric_factor/gmd_spectra); Spectra kernel removed from src/gmd_fast.cpp; RSpectra dropped from DESCRIPTION; eigencore 1.3.0 installed in default lib |
| 2 maxeig | done | gmdLA uses .metric_factor (Cholesky/eigen), no truncation; maxeig (default 5000) guards only the dense eigen of a singular general metric; warn_approx ignored |
| 3 remedy | done | default "error" in genpca/gpca_mle/PLS family; repair_metric() exported with report; genpca_metric_repaired warning class; version 0.2.0; vignettes/README/NEWS updated. rpls untouched (own penalty factorization). plsutils::partial_eig_once still pmax-clamps (dead helper, only tests call it) |
| 4 cov split | done | genpca_cov() GMD-only; new exported geigen_cov(); method = "geigen" deprecated alias with warning; tests re-pointed |
| 5 gpca_mle | done | default scale_fix = "none"; loglik recomputed at returned metrics; loglik_unpenalized and loglik_rescale_delta added; docs rewritten; hardening test updated; audit 8.7 active |
| 6 sparse | docs done; code partial | gpca-scale vignette documents what stays sparse; get_chol_lower_dense() refuses sparse input; spectra path already uses sparse CHOLMOD factors and never factors a singular large-side metric. Deferred: a user-facing pre-factored-metric input class (would need %*% methods everywhere the metric is multiplied) |
| 7 docs | done | NEWS 0.2.0 for phases 1-5 + swap; README, vignettes, WORDLIST updated; R CMD check --as-cran: 0 errors, 0 notes, 1 pre-existing machine-local clang-flag warning (R's Boolean.h); spelling: only pre-existing British spellings in untouched vignettes remain |

Review of Phases 1-2 + swap (fresh-context, 2026-09-04): numerics confirmed
across ~150 cases; 1 blocker fixed (gmdLA dual-branch norm guard scaled with
the data; now dimensionless), should-fixes applied (spectra no longer factors
a singular large-side metric: symmetric small-side form; sparse Cholesky pivot
check; explicit solver tol 1e-10; RSpectra skips removed from tests;
warn_approx deprecation warning; auto routes singular small-side metric above
maxeig to deflation; is_spd default tol = sqrt(eps); genpls alias mapping).
Added tests: 8.4 dual/up-scaling, test-solver-backend.R. Not done: uniform
component filter in genpca() (each backend filters with rank_rtol; equivalent).

Review of Phases 3-5 (fresh-context, 2026-09-04): Phases 3-4 verified end to
end (every remedy x storage type, warning class and report contents,
PLS-family threading, genpca_cov/geigen_cov split, vignettes all render).
1 blocker fixed: gpca_mle's covariance updates went through the tolerant
ensure_spd(), whose relative PD margin (from Phase 1) fired once
max(diag(Sigma)) > 1e6 * lambda and replaced the learned metric by a
near-identity (non-monotone path, wrong M/A). Now .mle_symmetric_pd():
average triangles, repair only when not PD at rtol = 0. Regression test added
(seed 2, 30x6 and 100x scale). Should-fixes applied: geigen_cov() uses the
shared .prep_one_metric() (its private ensure_psd() repaired silently);
geigen_cov own @param C; lambda = 0 warning; genpca_cov() warns when a
non-error constraints_remedy is supplied (ignored by the GMD form); Rd item
braces; ignore_attr in tests. Noted, pre-existing: truncate.Rd codoc mismatch,
genpls.Rd non-ASCII dash.

Independent release-review remediation (2026-09-04): all six findings in
`2026-09-04-release-review.md` addressed. Positive diagonal metric support is
consistent; deflation has a separate rank tolerance and the common output
filter keeps component arrays aligned; raw zero-component C++ vectors are
materialized correctly; eigencore nonconvergence is rejected centrally;
explicit clip bypasses tolerant acceptance; projected geigen semantics and
manual equations are corrected. Added 101 regression expectations, regenerated
Rd/Rcpp bindings, excluded diagnostic build artifacts, removed global vignette
title-option mutations and replaced stale submission claims with actual checks.

Final exact-tarball checks with eigencore 1.0.3 and 1.3.0: 0 errors, 0 warnings,
2 notes (new submission; local old HTML Tidy). Both pass all 1,403 expectations,
examples, vignette rebuilds and the PDF manual. Check-local Apple clang 15 plus
the existing Fortran runtime; user compiler settings and installed libraries
unchanged. Receipts and tarball SHA in `release-review-2026-09-04/fixed/README.md`.
Linux/Windows/external-builder checks remain outstanding. No commits or pushes.
