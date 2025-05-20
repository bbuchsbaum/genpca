library(testthat)
library(Matrix)

context("prep_constraints and construct_spatial_penalty utils")

# Test 1: prep_constraints error on non-PSD matrix with remedy="error"
test_that("prep_constraints errors for non-PSD constraints when remedy=error", {
  set.seed(1)
  X <- matrix(rnorm(4), 2, 2)
  A_bad <- Matrix(c(-0.5, 0, 0, 1), 2, 2, sparse=TRUE)
  expect_error(prep_constraints(X, A_bad, NULL, remedy="error"),
               "positive semi-definite")
})

# Test 2: prep_constraints with remedy="ridge" fixes negative eigenvalues

test_that("prep_constraints ridge remedy produces PSD matrix", {
  set.seed(1)
  X <- matrix(rnorm(4), 2, 2)
  A_bad <- Matrix(c(-0.5, 0, 0, 1), 2, 2, sparse=TRUE)
  pc <- prep_constraints(X, A_bad, NULL, remedy="ridge")
  eig <- eigen(as.matrix(pc$A), symmetric=TRUE, only.values=TRUE)$values
  expect_true(min(eig) >= -1e-8)
  expect_equal(dim(pc$A), c(2,2))
})

# Test 3: construct_spatial_penalty edge cases

test_that("construct_spatial_penalty handles knn edge cases", {
  spat_cds <- matrix(runif(2*5), nrow=2)
  expect_warning(Omega <- construct_spatial_penalty(spat_cds, k=5),
                 "Reducing k")
  expect_equal(dim(Omega), c(5,5))
  expect_true(Matrix::isSymmetric(Omega))
  expect_error(construct_spatial_penalty(spat_cds, k=0), "at least 1")
})

