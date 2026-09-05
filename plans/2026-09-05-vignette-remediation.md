# Vignette remediation — 2026-09-05

Implemented the reader/correctness audit across all six vignettes and the
pkgdown article index. No package implementation changes.

- Getting started now leads with inverse-variance-weighted USArrests, exposes
  scores/components/reconstruction/projection, and extracts all available
  components for the scree example. Theory and historical context follow the
  practical workflow. Corrected covariance/precision and row/column orientation
  guidance, preprocessing notation, singular-metric qualifications, and labels.
- Metric recipes distinguish alternating marginal scaling from covariance
  estimation. The experimental MLE example uses the default scale convention
  and shows eigenvalues, condition diagnostics, and likelihood progress rather
  than claiming identity recovery from a heatmap. Corrected repair and scalar
  normalization explanations.
- Scale guidance describes current factorization, dense data copies, metric
  guards, fill-in, and fallback costs. The randomized comparison reports its
  observed relative error instead of claiming plotting-precision agreement.
- Structured-noise conclusions are limited to the demonstrated simulations;
  removed universal separation and worst-case claims, and clarified that
  learning a metric pair still yields a shared pair.
- Sparse PCA plots align displayed signs to truth and share a color scale with
  an explicit zero. Fitted factors remain unchanged for reconstruction. Support
  counts acknowledge imperfect recovery.
- PLS quickstart uses related blocks and projects both into latent coordinates;
  the dense reference remains a separate contributor section.
- Added direct article links, ordered navigation, visible fitting warnings,
  and quiet checks on key example dimensions and recovery claims.

## Verification

- All six current-source vignettes rendered in separate R processes with
  `pkgload::load_all(export_all = FALSE)`, `LC_ALL=C`, and `LANG=C`.
- Revised figures visually inspected; numerical checks on dimensions,
  finiteness, projection agreement, signal recovery, and metric spectra passed.
- Matching title/index metadata, internal HTML links/anchors, and absence of
  visible unexpected warning output verified. The deliberate sfpca input-shape
  error remains part of the teaching example.
- `R CMD build` installed the package and rebuilt all six vignettes successfully.
  The tarball contains byte-identical current Rmd sources and six built HTMLs.
- The pkgdown article index built successfully to a scratch destination, using
  a scratch R cache; verified the three intended article groups.
- `git diff --check` passed. Full package tests and `--as-cran` were not repeated
  for this documentation-only change; this receipt does not replace the prior
  release check record.

Artifacts: `/tmp/genpca-vignette-fixed/` (renders, objects, figure extracts,
render logs, build log, verification script, and local article index).

Built package: `/tmp/genpca-vignette-fixed/genpca_0.2.0.tar.gz`

SHA-256: `2946b790b1a26da803ce0df74913fe3736e61ef66b7e480ac31f0f321aae79fa`

No site publication, commit, or push performed in this step.
