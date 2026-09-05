# Independent CRAN release review, 2026-09-04

Remediation follow-up: all six findings below have been fixed; see
[final checks and receipts](release-review-2026-09-04/fixed/README.md).
The original review below records the pre-fix candidate.

**Decision: blocked; do not submit or commit this as release-ready.** The
remediation improves the package substantially, but the final working tree has
two reproducible numerical defects, incomplete convergence handling, and broken
manual equations. No implementation files were changed or committed during this
review. Reproductions and check receipts are in `release-review-2026-09-04/`.

## Findings

### R1 — P1: positive singular values can have zero loading vectors

Locations: `R/gmd_fast.R:49-56`, `R/gmd_fast.R:143-146`,
`R/gpca_cov.R:198-205`; the same diagonal forward/inverse mismatch appears in
`R/internal_ops.R:85-86` and `R/weight_operators.R:29-34`.

The diagonal factors retain every positive weight in the forward square root,
but set the inverse square root to zero below `metric_rtol * max(weight)`.
Thus the target includes directions that cannot be mapped back to loadings.

With `X <- diag(c(1, 1e6))`, `A <- diag(c(1, 1e-10))`, `ncomp = 1`, and default
preprocessing, both `genpca(method = "eigen")` and `genpca(method = "spectra")`
return `sdev = 10` and `ov = c(0, 0)`. `genpca_cov(crossprod(X), R = A,
ncomp = 1)` likewise returns `d = 10`, `v = c(0, 0)` and `R_rank = 1`.
The defining invariant `t(ov) %*% A %*% ov = 1` instead equals zero. The eigen
fit also has zero scores. Deflation and randomized return a nonzero normalized
loading for this example. This reproduces under eigencore 1.0.3 and 1.3.0.

Fix: choose one numerical support consistently for all forward and inverse
operators and component extraction. If small positive weights are retained in
the objective, retain their inverses; if they are deliberately discarded,
remove those directions before solving and make that numerical policy explicit.
Test metric orthonormality, reconstruction and backend agreement on data with
large variance in a small-weight direction, including row weights and PLS.

Evidence: `numerical-probes.R` and both `numerical-probes-eigencore-*.log` files.

### R2 — P1: deflation ignores the new public rank_rtol

Locations: `R/gpca.R:451-495`, `R/gpca.R:1217`, `R/gpca.R:1298-1301`,
`src/gpca.cpp:83,137`.

The deflation calls pass only `threshold` as `thr`; neither implementation
receives `rank_rtol`, and the common return path does not filter it. For
`X <- diag(c(10, 1, .1))`, `ncomp = 3`, `rank_rtol = .2`, both C++ and R
deflation return all three singular values. Eigen, spectra and randomized
return only 10. This contradicts the public help and NEWS claim that all
methods apply `d_j > rank_rtol * d_1`. Auto selection can expose the same issue.

Fix: separate iteration tolerance from component acceptance in both deflation
kernels. Apply the same final component filter to every backend and keep all
vectors, variance summaries and component counts aligned. Simply filtering
afterwards is insufficient when a smaller `rank_rtol` asks to retain a
component that the existing `thr` stopping rule already discarded.

Evidence: `numerical-probes.R`, reproduced under both eigencore versions.

### R3 — P1: new equations fail manual rendering

Locations: `R/gpca.R:199`, `R/gpca_cov.R:33`, and generated
`man/genpca.Rd:105`, `man/genpca_cov.Rd:41`, `man/geigen_cov.Rd:36`.

The new equations nest `\code{rank_rtol}` inside `\eqn{...}`. Full
`R CMD check --as-cran` reports `Missing $ inserted` and related LaTeX errors
at these equations, plus HTML math errors: `Expected 'EOF', got '_'`.
The problem is present with both tested eigencore versions.

Fix the roxygen equations using valid math notation, or put the comparison in
plain code outside `\eqn`, regenerate Rd, and validate both PDF and HTML.
The fallback manual build additionally encounters the local missing font
`pcrr8t`; that environment problem must be distinguished from the independently
reported package equation errors.

Evidence: `check-eigencore-*.log`; scratch fallback log at
`/tmp/genpca-release-review/genpca.Rcheck/Rdlatex.log`.

### R4 — P2: convergence is reported by wrappers but ignored by several callers

Locations: `R/gplssvd_op.R:173-178`, `R/gpca.R:1024-1025`,
`R/gpca_cov.R:229`, `R/mnpca_mrl.R:633`.

`.top_svd()` and `.top_eigs_sym()` return `nconv` and `converged`, but these
callers use the values/vectors without checking them. For the seeded 110 x 80
and 110 x 70 matrices in `convergence-probe.R`, `tol = 1e-16` produces
`nconv = 0`, `converged = FALSE` under both eigencore versions.
`gplssvd_op(..., svd_opts = list(tol = 1e-16))` nevertheless returns an ordinary
result with no warning or convergence field. Values happen to be accurate in
this example; it proves discarded solver status, not large numerical error.
The spectra GMD paths and SFPCA already check status correctly.

Fix: consistently reject, warn, or invoke a bounded fallback for unconverged
results. Cover incomplete-result handling, not just successful convergence.

Evidence: `convergence-probe.R` and `convergence-eigencore-*.log`.

### R5 — P2: public explicit clip bypasses the strict clipping implementation

Locations: `R/repair_metric.R:65-73`, `R/gpca.R:71-80`, `NEWS.md:25-26`.

`repair_metric(diag(c(1, -1e-10)), method = "clip")` returns the negative
eigenvalue unchanged and reports `changed = FALSE`. The tolerant `is_psd`
gate prevents `clip_psd()` from running. General near-PSD metrics also bypass
the requested clip through `.prep_one_metric()`. Thus the internal strict
implementation does not establish the public NEWS guarantee that clip returns
a PSD result. The public `rtol` description explicitly permits returning the
input, creating a contradictory contract.

Fix: make an explicit clip request invoke strict clipping even when the input
passes tolerant validation, and report the actual change. Otherwise withdraw
the strict public guarantee and document the weaker semantics consistently.

Evidence: `public-contract-probes.R` and its log.

### R6 — P2: new geigen_cov documentation overstates its singular-metric equation

Locations: `R/gpca_cov.R:120-127`, `R/gpca_cov.R:323-341`.

The implementation solves a projected problem in `range(R)`, not necessarily
`C v = lambda R v` as advertised. With `C <- matrix(c(2,1,1,2), 2)` and
`R <- diag(c(1,0))`, it returns `v = c(1,0)`, `lambda = 2`; the full-equation
residual is `c(0,1)`. The projection is a defensible estimator, but its equation
and additional range constraint need to be explicit. The underlying range
restriction predates this change; exporting it under a new standalone name
is the point to correct its public contract.

Also clarify the comparison with GMD: commutation does not ensure the same
leading component or ordering. For `C = diag(c(4,1))`, `R = diag(c(9,1))`, GMD
selects the first coordinate and geigen selects the second.

Fix: document the projected equation and range constraint, or require the
compatibility conditions needed for the full equation. Add an independent
equation-residual test and a diagonal example where component order differs.

Evidence: `public-contract-probes.R` and its log.

## Release checks and provenance

- Package: genpca 0.2.0; reviewed dirty tree based on
  `1a0ff8de9108bb9c8cd66ff3cb98fd41c39ed004`.
- Built source: `/tmp/genpca_0.2.0.tar.gz`, SHA-256
  `4c737e9c643299a801b53585c2bb54e52b0155b8e3c5b31c78a986ac4c1ad895`.
- `R CMD build` rebuilt vignettes successfully. Both full checks used this
  exact tarball with `LC_ALL=C LANG=C`, without `--no-manual`.
- Platform: macOS 14.3 arm64, R 4.5.1, Matrix 1.7-3; compilation used Homebrew
  clang 20.1.8. Installed multivarious was 0.3.1. This is not an all-CRAN or
  cross-platform dependency matrix.
- Full check with installed eigencore 1.3.0: **1 error, 2 warnings, 4 notes**.
- Full check with CRAN eigencore 1.0.3 in an isolated scratch library:
  **1 error, 2 warnings, 4 notes**. The user's installed library was unchanged.
- Each check passed tests: **1,302 expectations, 0 failures, 74 test warnings,
  1 empty-test skip**. These testthat warnings are distinct from check-level
  warnings. Examples and vignette rebuilds passed with both versions.
- Check warnings: local R-header clang flag warning, PDF manual equations.
  Error: fallback PDF manual generation. Notes: new submission, stray
  `Rplots.pdf`, HTML validation/math diagnostics, leftover `genpca-manual.tex`
  from manual failure. The HTML note also identifies the local old `tidy`.
- `git diff --check` passed before review-only artifacts were added.
- Live CRAN package index retrieved with `available.packages()` listed
  eigencore 1.0.3, multivarious 0.3.2, albersdown 2.0.0 and adjoin 0.1.0;
  genpca was absent. The eigencore minimum dependency is available on CRAN and
  passed the existing suite here. Its development version is not required to
  reproduce these findings.

## Packaging and submission follow-up

The built tarball actually contains `vignettes/Rplots.pdf`,
`tests/testthat/_problems/` and `tests/testthat/testthat-problems.rds`. Exclude or
remove these exact artifacts before rebuilding; this review preserved them.
`.DS_Store` was untracked but not present in the tarball.

`cran-comments.md` still names RSpectra headers and lists win-builder and
macOS builder environments without fresh receipts for this dirty candidate.
Update it from the final candidate's evidence. The repository currently has
a pkgdown workflow but no R CMD check workflow. Obtain Linux and Windows
release/devel validation, including appropriate minimum-version coverage,
before submission. Neither local checks nor older hosted results establish
those gates for this candidate.

The deferred prefactored-metric interface and upstream eigencore performance
work are separate from the blockers above. The sparse limitations are largely
documented honestly: the public spectra and randomized routes still densify X.

Only this report and its evidence directory were added during review. No fixes,
deletions of the pre-existing artifacts, commit, push, or submission were made.
