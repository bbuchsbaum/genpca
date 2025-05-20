context("solve_gep_subspace algorithmic correctness")

library(testthat)
library(Matrix)

# Helper to check generalized eigen equation: S1 V = S2 V diag(lambda)
check_gep <- function(S1, S2, V, lambda, tol=1e-6) {
  left <- S1 %*% V
  right <- S2 %*% (V %*% diag(lambda))
  expect_equal(left, right, tolerance = tol)
  gram <- t(V) %*% S2 %*% V
  expect_equal(gram, diag(length(lambda)), tolerance = tol)
}

set.seed(123)

# Small dimension for exact comparison
n <- 8
q <- 3

# Generate random SPD matrices
A <- matrix(rnorm(n * n), n, n)
S2 <- crossprod(A) + diag(n) * 0.1
B <- matrix(rnorm(n * n), n, n)
S1 <- crossprod(B) + diag(n) * 0.1

# Cast to Matrix class for consistency
S1 <- Matrix(S1)
S2 <- Matrix(S2)

# Direct eigen decomposition of solve(S2) %*% S1
M <- solve(S2, S1)
E <- eigen(M)

true_largest_vals <- E$values[1:q]

# Test largest eigenvalues
res_largest <- genpca:::solve_gep_subspace(S1, S2, q = q, which = "largest",
                                           reg_S = 1e-6, reg_T = 1e-8)

test_that("solve_gep_subspace matches direct solution for largest eigenvalues", {
  expect_equal(res_largest$values, true_largest_vals, tolerance = 1e-5)
  check_gep(S1, S2, res_largest$vectors, res_largest$values, tol = 1e-5)
})

# For smallest eigenvalues
true_smallest_vals <- rev(E$values)[1:q]
res_smallest <- genpca:::solve_gep_subspace(S1, S2, q = q, which = "smallest",
                                            reg_S = 1e-6, reg_T = 1e-8)

test_that("solve_gep_subspace matches direct solution for smallest eigenvalues", {
  expect_equal(res_smallest$values, true_smallest_vals, tolerance = 1e-5)
  check_gep(S1, S2, res_smallest$vectors, res_smallest$values, tol = 1e-5)
})

