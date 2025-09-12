testthat::test_that("deflation closely matches eigen on small dense problems", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(123)
  n <- 60; p <- 40; k <- 8
  X <- matrix(rnorm(n * p), n, p)

  # Eigen (reference) vs Deflation (R path)
  fit_eig <- genpca(X, ncomp = k, method = "eigen",
                    preproc = multivarious::center())
  fit_def <- genpca(X, ncomp = k, method = "deflation", use_cpp = FALSE,
                    preproc = multivarious::center(), threshold = 1e-8)

  # 1) Singular values: tight numerical match
  testthat::expect_equal(fit_eig$sdev[1:k], fit_def$sdev[1:k], tolerance = 1e-4)

  # 2) Subspace agreement via principal angles (scores space)
  # Use ou (orthonormal in M metric). With identity metrics here, it is Euclidean-orthonormal.
  Ou_e <- as.matrix(fit_eig$ou)[, 1:k, drop = FALSE]
  Ou_d <- as.matrix(fit_def$ou)[, 1:k, drop = FALSE]
  sv  <- svd(t(Ou_e) %*% Ou_d)$d
  # Sum of squared canonical correlations should be ~ k
  testthat::expect_equal(sum(sv^2), k, tolerance = 1e-3)

  # 3) Scores and loadings: allow arbitrary orthogonal rotation within the subspace.
  procrustes_align <- function(B, A) {
    # Find orthogonal R minimizing ||A - B R||_F via SVD(B^T A).
    BtA <- crossprod(B, A)
    sv  <- svd(BtA)
    R   <- sv$u %*% t(sv$v)
    B %*% R
  }

  S_e <- as.matrix(multivarious::scores(fit_eig))[, 1:k, drop = FALSE]
  S_d <- as.matrix(multivarious::scores(fit_def))[, 1:k, drop = FALSE]
  S_d_al <- procrustes_align(S_d, S_e)
  testthat::expect_lt(norm(S_e - S_d_al, "F"), 1e-3 * (1 + norm(S_e, "F")))

  L_e <- as.matrix(multivarious::components(fit_eig))[, 1:k, drop = FALSE]
  L_d <- as.matrix(multivarious::components(fit_def))[, 1:k, drop = FALSE]
  L_d_al <- procrustes_align(L_d, L_e)
  testthat::expect_lt(norm(L_e - L_d_al, "F"), 1e-3 * (1 + norm(L_e, "F")))
})
