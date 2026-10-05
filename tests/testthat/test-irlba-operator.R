test_that("matrix-free irlba products preserve orientation and block shapes", {
  for (dims in list(c(9L, 7L), c(7L, 9L))) {
    A <- matrix(seq_len(prod(dims)) / 17, dims[1], dims[2])
    calls <- c(forward = 0L, adjoint = 0L)
    op <- .irlba_operator(list(dims = dims,
      S_mv = function(v) {
        calls["forward"] <<- calls["forward"] + 1L
        A %*% v
      },
      ST_mv = function(u) {
        calls["adjoint"] <<- calls["adjoint"] + 1L
        crossprod(A, u)
      }))
    v <- seq_len(ncol(A)) / 5
    u <- seq_len(nrow(A)) / 3
    V <- cbind(v, -v)
    U <- rbind(u, -u)
    expect_identical(dim(op), dims)
    expect_equal(op %*% v, A %*% v)
    expect_equal(u %*% op, u %*% A)
    expect_equal(op %*% V, A %*% V)
    expect_equal(U %*% op, U %*% A)
    expect_equal(drop(u %*% (op %*% v)), drop((u %*% op) %*% v))
    expect_true(all(calls > 0L))
    expect_error(op %*% numeric(ncol(A) + 1L), "non-conformable")
    expect_error(numeric(nrow(A) + 1L) %*% op, "non-conformable")
  }
})

test_that("irlba evaluates both matrix-free products on rectangular operators", {
  skip_if_not_installed("irlba")
  for (dims in list(c(80L, 70L), c(70L, 80L))) {
    # Deterministic nonzero operator with separated leading singular values.
    A <- matrix(sin(seq_len(prod(dims))), dims[1], dims[2])
    diag(A) <- diag(A) + seq_len(min(dims)) / 2
    calls <- c(forward = 0L, adjoint = 0L)
    op <- .irlba_operator(list(dims = dims,
      S_mv = function(v) {
        calls["forward"] <<- calls["forward"] + 1L
        A %*% v
      },
      ST_mv = function(u) {
        calls["adjoint"] <<- calls["adjoint"] + 1L
        crossprod(A, u)
      }))
    fit <- irlba::irlba(op, nv = 3, nu = 3, tol = 1e-9)
    expect_true(all(calls > 0L))
    expect_equal(fit$d, svd(A, nu = 0, nv = 0)$d[1:3], tolerance = 1e-7)
    expect_identical(dim(fit$u), c(dims[1], 3L))
    expect_identical(dim(fit$v), c(dims[2], 3L))
    expect_equal(A %*% fit$v, sweep(fit$u, 2, fit$d, `*`), tolerance = 1e-7)
    expect_equal(crossprod(A, fit$u), sweep(fit$v, 2, fit$d, `*`), tolerance = 1e-7)
  }
})


test_that("skinny irlba operators cannot silently enter a dense fallback", {
  skip_if_not_installed("irlba")
  X <- matrix(sin(seq_len(20L * 80L)), 20L, 80L)
  Y <- matrix(cos(seq_len(20L * 5L)), 20L, 5L)
  expect_error(gplssvd_op(X, Y, k = 2, svd_backend = "irlba"),
               "requires both column counts >= 6")
  expect_error(gplssvd_op(Y, X, k = 2, svd_backend = "irlba"),
               "requires both column counts >= 6")
})
