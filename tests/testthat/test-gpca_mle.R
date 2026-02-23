library(testthat)
library(genpca)

test_that("gpca_mle returns SPD metrics with correct dimensions", {
  set.seed(1)
  X <- matrix(rnorm(60), 10, 6)
  res <- suppressWarnings(
    gpca_mle(X, ncomp = 3, max_iter = 4, lambda = 1e-3, verbose = FALSE)
  )
  expect_equal(dim(res$A), c(6, 6))
  expect_equal(dim(res$M), c(10, 10))
  expect_true(genpca:::is_spd(res$A))
  expect_true(genpca:::is_spd(res$M))
  expect_s3_class(res$fit, "genpca")
  expect_true(length(res$loglik_path) >= 1)
  expect_true(is.finite(res$loglik))
})

test_that("gpca_mle produces a finite likelihood path", {
  set.seed(2)
  X <- matrix(rnorm(40), 8, 5)
  res <- suppressWarnings(
    gpca_mle(X, ncomp = 2, max_iter = 5, lambda = 1e-3, verbose = FALSE)
  )
  expect_gt(length(res$loglik_path), 1)
  expect_true(all(is.finite(res$loglik_path)))
})
