library(testthat)

# Tests for explicit column-likeness and automatic path selection


test_that("force_col_likeness chooses column path and matches direct SVD", {
  set.seed(1)
  n <- 8; p <- 5; q <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  fit_direct <- genplsc(X, Y, ncomp = 2, dual = FALSE, verbose = FALSE)
  fit_col    <- genplsc(X, Y, ncomp = 2, dual = TRUE,
                        force_col_likeness = TRUE, verbose = FALSE)

  expect_false(fit_col$want_row_like)
  expect_equal(fit_direct$d_full, fit_col$d_full, tolerance = 1e-2)

  cor_vx <- abs(cor(as.vector(fit_direct$vx), as.vector(fit_col$vx)))
  cor_vy <- abs(cor(as.vector(fit_direct$vy), as.vector(fit_col$vy)))
  expect_gt(cor_vx, 0.9)
  expect_gt(cor_vy, 0.9)
})


test_that("automatic row/column likeness decision follows n < (p+q)/2 rule", {
  set.seed(2)
  n_small <- 3; p <- 4; q <- 5  # n_small < (p+q)/2 => row likeness
  Xs <- matrix(rnorm(n_small * p), n_small, p)
  Ys <- matrix(rnorm(n_small * q), n_small, q)
  fit_small <- genplsc(Xs, Ys, ncomp = 1, dual = TRUE, verbose = FALSE)
  expect_true(fit_small$want_row_like)

  n_big <- 20  # n_big > (p+q)/2 => column likeness
  Xb <- matrix(rnorm(n_big * p), n_big, p)
  Yb <- matrix(rnorm(n_big * q), n_big, q)
  fit_big <- genplsc(Xb, Yb, ncomp = 1, dual = TRUE, verbose = FALSE)
  expect_false(fit_big$want_row_like)
})


test_that("force_row_likeness overrides automatic column choice", {
  set.seed(3)
  n <- 15; p <- 4; q <- 5  # this would normally trigger column likeness
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  fit_auto <- genplsc(X, Y, ncomp = 1, dual = TRUE, verbose = FALSE)
  expect_false(fit_auto$want_row_like)

  fit_forced <- genplsc(X, Y, ncomp = 1, dual = TRUE,
                        force_row_likeness = TRUE, verbose = FALSE)
  expect_true(fit_forced$want_row_like)
})

