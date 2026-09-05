# Replacing RSpectra with eigencore in genpca: findings (2026-09-04)

Status: evaluation complete for the R-level call sites; recommendation and
plan below. Nothing in the package has been changed. All scripts and the
experimental shim are in `plans/eigencore-swap/` (build-ignored).

Versions tested: eigencore 1.0.3 (CRAN, installed to a scratch lib) and
eigencore 1.3.0 (dev source at ~/code/eigencore, installed to a scratch lib);
RSpectra 0.16.2. The copy of the library in the default R library is a stale
1.0.0.9000 and should be reinstalled.

## 1. Where genpca uses RSpectra

R level (12 call sites, all covered by the shim):

| Site | Call | Notes |
|---|---|---|
| `gmd_fast.R:39` | `svds(Xw, k, opts = list(maxitr, tol))` | diagonal-metric spectra path |
| `plsutils.R:44` | `eigs_sym(M, k, which, opts = list(tol))` | |
| `sfpca.R:541` | `svds(A = fn, k = 1, Atrans = fn, dim = c(n, p))` | function operator, has base-svd fallback |
| `sfpca.R:565` | `eigs_sym(dgC, k = 1, "LM")` | |
| `mnpca_mrl.R:633` | `svds(Y, k, nu, nv)` | |
| `gplssvd_op.R:168` | `svds(A = fn, k, nu, nv, opts, Atrans, dim)` | function operator, no fallback |
| `gpca_cov.R:219` | `eigs_sym(B, k, "LA")` | |
| `gpca_cov.R:280,323` | `eigs_sym(C, 1, "SA")` | min eigenvalue, p > 800 |
| `gpca.R:959` | `eigs_sym(sparse metric, k = maxeig, "LM")` | the 8.1 truncation branch; deleted in Phase 2 of the audit plan |
| `gpca.R:1001` | `eigs_sym(target, k, "LM")` | main eigen-path top-k |
| `tests/.../test_sfpca_cd.R:27` | direct `RSpectra::eigs_sym` | test helper |

C++ level: `src/gmd_fast.cpp` includes `<Spectra/SymEigsSolver.h>` via
`LinkingTo: RSpectra` and runs `spectra_topk()` on matrix-free operators
`L_Rᵀ XᵀQX L_R` (primal) and `L_Qᵀ X R Xᵀ L_Q` (dual). eigencore has no
`inst/include` and no `R_RegisterCCallable` entry points, so this cannot be
swapped at the C++ level. See section 4 for the way around it.

## 2. Equivalence results

Accuracy (both versions identical on these):

| Case | eigencore vs LAPACK | RSpectra vs LAPACK |
|---|---|---|
| dense sym 300, k=5, LM values | 4e-15 | 2e-14 |
| dense sym 300, k=5, LA values | 2e-15 | 2e-14 |
| k=1 SA (min eigenvalue, indefinite) | exact match | exact match |
| subspace angle, k=5 (default tol) | 2.4e-7 | 2.6e-8 |
| svds dense 400x60, k=5, d | 2e-14 | 6e-14 |
| svds u/v subspace | 3.7e-8 | 1.5e-8 |
| matrix classes: base, dgeMatrix, dsyMatrix, dgCMatrix, dsCMatrix | all accepted | |
| k = n, k = min(n,p) | dense fallback, exact | warning / error |
| clustered eigenvalues split by k | exact, certificate passed | |

Full genpca test suite with every RSpectra call routed through the shim:

| Library | tests | failed | errors |
|---|---|---|---|
| master baseline (RSpectra) | 183 | 2 | 0 |
| eigencore 1.0.3 | 183 | 3 | 0 |
| eigencore 1.3.0 | 183 | 3 | 0 |

The 2 baseline failures are pre-existing rpls message-regex expectations in
`test-coverage-95-branches.R:558,570`, unrelated to the solver. The one extra
failure is an artifact of the experiment, not of eigencore:
`test-coverage-gpca-branches.R:159` matches the literal warning text
"using RSpectra", and the blanket text substitution that built the scratch
copy rewrote that string inside the warning message too. With a real wrapper
the warning text is untouched and the suite result equals the baseline. (The
branch it tests is deleted in Phase 2 of the audit plan anyway.) Every numerical equivalence test
(eigen vs spectra, deflation vs eigen, CCA vs genpls, genpls vs gplssvd_op,
dual-path, randomized) passes unchanged at its existing tolerance.

GMD formulated as an SVD of the whitened operator (section 4), 4000x600,
general dense SPD column metric, diagonal row metric, k = 10, against the
C++ Spectra kernel (`genpca(method = "spectra")`): sdev relative difference
5.7e-15, V subspace angle 2.9e-7. Timing in section 4.

## 3. Behavioural differences a wrapper must handle

1. `opts` is accepted and ignored. Map `opts$tol` → `tol` and `opts$maxitr`
   → `maxit` (eigen) yourself.
2. No `Atrans`/`dim`/`args` for function operators. Wrap with
   `linear_operator(dim, apply, apply_adjoint)`.
3. Callbacks are block GEMM-style: `apply(X, alpha = 1, beta = 0, Y = NULL)`
   with `X` an n x b matrix (b = 1 or 2 observed) and must return a double
   matrix. Observed alpha, beta are always (1, 0). A callback returning a
   numeric vector, or one whose first formal is not named `X` on 1.0.3
   (1.0.3 passes it by name, 1.3.0 positionally), fails with the opaque
   "native matrix-free Golub-Kahan failed with status=-8".
4. `.Random.seed` is consumed on every call, even with `seed =`. Values are
   bitwise reproducible across calls; vectors reproduce only to ~3e-8 (same
   with or without a fixed seed, so `seed` does not pin the start on the
   native paths). The shim saves and restores `.Random.seed`.
5. Certificates: `passed` is withheld (FALSE) for every matrix-free operator
   because the norm bound is a Hutchinson estimate. Check `nconv`,
   `converged` and residuals for operator paths rather than `passed`.
6. Matrix-free symmetric eigenproblems (`eigs_sym` on a `linear_operator`
   with `hermitian()`): the default `auto()` plan falls to an R-level
   "prototype/oracle" Lanczos in both versions; it returns values wrong at
   1e-5 to 1e-2 with `nconv = 0` and no R-level warning. Only 1.3.0 with
   `method = lanczos(block >= 2)` runs the native block path (accurate to
   1e-14). genpca should not use this path at all; formulate GMD as an SVD.
7. Default `tol` is 1e-8 (RSpectra 1e-10). Pass `tol` explicitly where
   genpca's tests assume tight subspaces.

## 4. Performance

Clean sequential runs, median of 3, nothing else running, installed
libraries (`plans/eigencore-swap/timing.R`). Planner method names are what
eigencore printed.

| Problem | RSpectra | eigencore 1.0.3 | eigencore 1.3.0 | eigencore method |
|---|---|---|---|---|
| dense sym 1500, k=10 | 0.038 s | 0.078 s (0.075 uncertified) | 0.366 s (0.346 uncertified) | native scalar thick-restart Hermitian Lanczos |
| sparse sym 20000, ~10 nnz/row, k=10 | 0.37 s | 0.40 s | 0.86 s | native scalar thick-restart Hermitian Lanczos |
| svds dense 5000x300, k=10 | 0.020 s | 0.015 s | 0.202 s | native certified Gram SVD special case |
| operator svds 400x300 (gplssvd-type), k=5 | 0.15 s | 0.32 s | 0.31 s | native matrix-free Golub-Kahan callback cycle |
| GMD operator 4000x600, general metric, k=10 | 0.43 s (C++ Spectra kernel) | 0.76 s | 0.65 s | native matrix-free Golub-Kahan callback cycle |
| per-call overhead at test sizes | ~0 | 2-4 ms | 2-4 ms | |

Reading: on 1.0.3 eigencore is at parity with RSpectra on sparse and dense
SVD, 2x slower on dense symmetric eigen, and about 2x slower on operator
problems because each block matvec crosses the R callback boundary. The
operator-formulated GMD is 1.5-1.75x slower than the in-process C++ Spectra
kernel; acceptable for a dogfood swap, and the gap is the callback boundary,
not the algorithm. **1.3.0 regresses 5-13x against 1.0.3 on the built-in
dense and sparse paths with identical planner method names**; fix that
upstream before pointing genpca at it.

## 5. Recommendation

Swap is feasible and numerically safe at the R level today with a thin
wrapper. Do it as part of the audit plan, not before it:

- **Phase 1a (new, after audit Phase 1).** Add `R/solver_backend.R` with
  `.top_eigs_sym(A, k, which, tol, maxit)` and
  `.top_svd(A, k, nu, nv, tol, adjoint, dim)` wrapping eigencore as in the
  shim (opts mapping, `linear_operator` for closures, GEMM-shaped callbacks,
  `.Random.seed` guard, explicit `tol`). Route the 11 R call sites through
  it. Keep the shim's `force()` (the closure captures `A` before it is
  rebound). Move `Imports: RSpectra` to `Suggests` only after step 1b.
- **Phase 1b.** Replace the C++ Spectra kernel in `gmd_fast.cpp` with an
  R-level `linear_operator` for the whitened matrix `B = L_Qᵀ X L_R`
  (apply: `v ↦ L_Qᵀ (X (L_R v))`, adjoint: `u ↦ L_Rᵀ (Xᵀ (L_Q u))`), solved
  with `svd_partial`; d = singular values, `ov = L_R^{-ᵀ} V`,
  `ou = L_Q^{-ᵀ} U`. This is exactly the GMD, avoids the weak symmetric
  matrix-free path, costs 1.5-1.75x the C++ kernel at 4000x600 (callback
  boundary), and lets the audit plan's Phase 6 sparse-factor work stay in R
  (`Matrix` triangular solves). `gmd_fast.cpp` then keeps only the dense-fallback and
  orthonormalization helpers; drop `LinkingTo: RSpectra, RcppEigen` if
  nothing else needs them.
- **Phase 1c.** Tests: keep every existing backend-equivalence test at its
  tolerance; add one test per item in section 3 (opts mapping, operator
  wrapping, RNG guard, `nconv` check), and a test that `svd_partial`'s
  certificate `converged` is TRUE on the GMD operator.
- **Upstream to eigencore (own repo):** (i) plan matrix-free hermitian
  problems onto the native block path or raise a warning when the prototype
  runs; (ii) the 1.3.0 regression: 5-13x slower than 1.0.3 on native dense
  Lanczos, sparse Lanczos and the Gram SVD special case; (iii) friendlier error
  than status -8 when a callback errors or returns a non-matrix; (iv) honour
  `seed` on native paths or document that it does not pin the start; (v)
  accept the first callback argument positionally in all versions (1.0.3
  requires the formal to be named `X`).
