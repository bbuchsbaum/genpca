test_that("gmd_fast_cpp dense/sparse dispatch works", {
  set.seed(42)
  n <- 100
  p <- 30
  k <- 5
  X <- matrix(rnorm(n * p), n, p)
  Qd <- diag(runif(n, 0.5, 2))
  Rd <- diag(runif(p, 0.5, 2))
  
  Qs <- methods::as(Matrix::Matrix(Qd, sparse = TRUE), "dgCMatrix")
  Rs <- methods::as(Matrix::Matrix(Rd, sparse = TRUE), "dgCMatrix")
  
  # Test dense version
  res_dense <- genpca:::gmd_fast_cpp(X, Qd, Rd, k)
  expect_equal(length(res_dense$d), k)
  expect_equal(nrow(res_dense$u), n)
  expect_equal(nrow(res_dense$v), p)
  expect_equal(ncol(res_dense$u), k)
  expect_equal(ncol(res_dense$v), k)
  
  # Test sparse version  
  res_sparse <- genpca:::gmd_fast_cpp(X, Qs, Rs, k)
  expect_equal(length(res_sparse$d), k)
  expect_equal(nrow(res_sparse$u), n)
  expect_equal(nrow(res_sparse$v), p)
  
  # Results should be very similar (same problem, different matrix types)
  expect_equal(res_dense$d, res_sparse$d, tolerance = 1e-6)
  
  # Check subspaces are aligned (may have sign/permutation differences)
  if (requireNamespace("clue", quietly = TRUE)) {
    align_u <- align_perm_sign(res_dense$u, res_sparse$u)
    expect_true(mean(align_u$corr) > 0.99)
    
    align_v <- align_perm_sign(res_dense$v, res_sparse$v)
    expect_true(mean(align_v$corr) > 0.99)
  }
})

test_that("gmd_fast_cpp handles dsyMatrix input", {
  set.seed(123)
  n <- 50
  p <- 20
  k <- 3
  X <- matrix(rnorm(n * p), n, p)
  
  # Create dense symmetric matrices (will become dsyMatrix)
  Q <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 0.1
  R <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p) * 0.1
  
  # Convert to Matrix format (will be dsyMatrix)
  Q_mat <- Matrix::Matrix(Q)
  R_mat <- Matrix::Matrix(R)
  
  expect_true(inherits(Q_mat, "dsyMatrix") || inherits(Q_mat, "dpoMatrix"))
  expect_true(inherits(R_mat, "dsyMatrix") || inherits(R_mat, "dpoMatrix"))
  
  # Should not error with dsyMatrix input
  res <- genpca:::gmd_fast_cpp(X, Q_mat, R_mat, k)
  
  expect_equal(length(res$d), k)
  expect_equal(nrow(res$u), n)
  expect_equal(nrow(res$v), p)
  expect_equal(res$k, k)
})