library(testthat)
library(Matrix)
library(RSpectra)

context("sfpca function")

test_that("Rank-1 matrix is recovered correctly", {
  set.seed(123)
  n <- 10
  p <- 5
  u_true <- rnorm(n)
  v_true <- rnorm(p)
  X <- u_true %*% t(v_true)
  spat_cds <- matrix(runif(p * 2), nrow = 2, ncol = p)
  result <- sfpca(X, K = 1, spat_cds = spat_cds, verbose = FALSE, lambda_v=.001, lambda_u=.001)
  
  u_est <- result$u[,1]
  v_est <- result$v[,1]
  d_est <- result$d[1]
  
  # Handle sign ambiguity
  sign_u <- sign(crossprod(u_est, u_true))
  u_est <- u_est * sign_u
  sign_v <- sign(crossprod(v_est, v_true))
  v_est <- v_est * sign_v
  
  # Check closeness
  expect_equal(u_est, u_true, tolerance = 1e-4)
  expect_equal(v_est, v_true, tolerance = 1e-4)
  expect_equal(d_est, crossprod(u_true, X %*% v_true), tolerance = 1e-4)
})

test_that("Orthogonal columns result in smooth components", {
  set.seed(123)
  n <- 100
  p <- 20
  X <- diag(p)[1:n,]
  spat_cds <- matrix(runif(p * 2), nrow = 2, ncol = p)
  result <- sfpca(X, K = 2, spat_cds = spat_cds, verbose = FALSE)
  
  v_est <- result$v
  
  # Compute smoothness measure based on spat_cds
  distances <- as.matrix(dist(t(spat_cds)))
  smoothness <- sum(distances * (v_est %*% t(v_est)))
  
  # Set a threshold for smoothness
  expect_true(smoothness < 100)  # Adjust threshold as needed
})

test_that("Sparse signals result in sparse components", {
  set.seed(123)
  n <- 50
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X[, sample(p, 10)] <- 0  # Introduce sparsity
  spat_cds <- matrix(runif(p * 3), nrow = 3, ncol = p)
  result <- sfpca(X, K = 2, spat_cds = spat_cds, lambda_u = 0.1, lambda_v = 0.1, verbose = FALSE)
  
  u_est <- result$u
  v_est <- result$v
  
  # Check sparsity
  expect_true(sum(abs(u_est) > 1e-4) < 20)  # Assuming sparsity threshold
  expect_true(sum(abs(v_est) > 1e-4) < 20)
})