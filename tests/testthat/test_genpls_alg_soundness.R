library(testthat)
library(Matrix)

# Algorithmic soundness tests for genpls (genplsr)

test_that("genpls loadings are orthonormal with identity constraints", {
  set.seed(100)
  n <- 8
  p <- 5
  q <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  fit <- genpls(X, Y, ncomp = 3, verbose = FALSE)

  K <- ncol(fit$vx)
  expect_equal(dim(fit$vx), c(p, K))
  expect_equal(dim(fit$vy), c(q, K))

  gram_x <- t(fit$vx) %*% fit$vx
  gram_y <- t(fit$vy) %*% fit$vy

  expect_equal(gram_x, diag(K), tolerance = 1e-6)
  expect_equal(gram_y, diag(K), tolerance = 1e-6)

  gram_tx <- t(fit$tilde_Tx) %*% fit$tilde_Tx
  gram_ty <- t(fit$tilde_Ty) %*% fit$tilde_Ty

  expect_equal(gram_tx, diag(K), tolerance = 1e-6)
  expect_equal(gram_ty, diag(K), tolerance = 1e-6)
})


test_that("genpls respects orthogonality under column constraints", {
  set.seed(200)
  n <- 10
  p <- 6
  q <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  Ax_tmp <- crossprod(matrix(rnorm(p * p), p, p))
  Ay_tmp <- crossprod(matrix(rnorm(q * q), q, q))
  Ax <- forceSymmetric(Ax_tmp)
  Ay <- forceSymmetric(Ay_tmp)

  fit <- genpls(X, Y, Ax = Ax, Ay = Ay, ncomp = 2, verbose = FALSE)

  K <- ncol(fit$vx)
  expect_equal(dim(fit$vx), c(p, K))
  expect_equal(dim(fit$vy), c(q, K))

  gram_x <- t(fit$vx) %*% (Ax %*% fit$vx)
  gram_y <- t(fit$vy) %*% (Ay %*% fit$vy)

  expect_equal(gram_x, diag(K), tolerance = 1e-6)
  expect_equal(gram_y, diag(K), tolerance = 1e-6)

  gram_px <- t(fit$tilde_Px) %*% fit$tilde_Px
  gram_py <- t(fit$tilde_Py) %*% fit$tilde_Py

  expect_equal(gram_px, diag(K), tolerance = 1e-6)
  expect_equal(gram_py, diag(K), tolerance = 1e-6)
})

