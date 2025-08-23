#' Fast generalized matrix decomposition (dense/sparse dispatch)
#' @param X numeric matrix (n x p)
#' @param Q,R constraints (weights/metrics) for rows/cols
#' @param k number of components
#' @param tol tolerance
#' @param maxit maximum iterations (for compatibility, ignored in new implementation)
#' @param seed random seed (for compatibility, ignored in new implementation)
#' @param topk logical; use top-k symmetric eigen when available (ARPACK). Defaults to TRUE.
#' @param cache logical; cache Cholesky factors across calls when dense. Defaults to TRUE.
#' @keywords internal
#' @importFrom methods as is
gmd_fast_cpp <- function(X, Q, R, k, tol = 1e-9, maxit = 1000L, seed = 1234L, topk = TRUE, cache = TRUE) {
  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)

  # coerce symmetric dense forms that trip Rcpp
  if (inherits(Q, "dsyMatrix")) Q <- methods::as(Q, "dgeMatrix")
  if (inherits(R, "dsyMatrix")) R <- methods::as(R, "dgeMatrix")

  n <- nrow(X); p <- ncol(X)
  primal <- (p <= n)  # use primal when small side is p

  # ---- choose path and apply caching for dense constraints ----
  if (primal) {
    if (!methods::is(R, "sparseMatrix") && isTRUE(cache)) {
      L_R <- get_chol_lower_dense(R)
      if (methods::is(Q, "sparseMatrix")) {
        res <- gmd_fast_cpp_primal_sp(X, methods::as(Q, "dgCMatrix"), L_R, k, tol, topk)
      } else {
        res <- gmd_fast_cpp_primal_dn(X, as.matrix(Q), L_R, k, tol, topk)
      }
    } else {
      # fall back to non-cached path
      if (methods::is(Q, "sparseMatrix") || methods::is(R, "sparseMatrix")) {
        res <- gmd_fast_cpp_sp(X, methods::as(Q, "dgCMatrix"), methods::as(R, "dgCMatrix"), k, tol, topk)
      } else {
        res <- gmd_fast_cpp_dn(X, as.matrix(Q), as.matrix(R), k, tol, topk)
      }
    }
  } else {
    # dual path (n < p): cache Cholesky for Q when dense
    if (!methods::is(Q, "sparseMatrix") && isTRUE(cache)) {
      L_Q <- get_chol_lower_dense(Q)
      if (methods::is(R, "sparseMatrix")) {
        res <- gmd_fast_cpp_dual_sp(X, L_Q, methods::as(R, "dgCMatrix"), k, tol, topk)
      } else {
        res <- gmd_fast_cpp_dual_dn(X, L_Q, as.matrix(R), k, tol, topk)
      }
    } else {
      if (methods::is(Q, "sparseMatrix") || methods::is(R, "sparseMatrix")) {
        res <- gmd_fast_cpp_sp(X, methods::as(Q, "dgCMatrix"), methods::as(R, "dgCMatrix"), k, tol, topk)
      } else {
        res <- gmd_fast_cpp_dn(X, as.matrix(Q), as.matrix(R), k, tol, topk)
      }
    }
  }

  # Ensure d is a vector (C++ functions may return it as a column matrix)
  if (is.matrix(res$d)) {
    res$d <- as.vector(res$d)
  }
  
  # Name outputs like multivarious - filter by tolerance first
  keep <- res$d > tol
  if (sum(keep) < length(res$d)) {
    res$u <- res$u[, keep, drop = FALSE]
    res$v <- res$v[, keep, drop = FALSE]
    res$d <- res$d[keep]
  }
  
  # IMPORTANT: The C++ implementation returns scores as U*D, but gpca.R expects orthonormal U
  # We need to normalize the scores to get orthonormal U
  if (length(res$d) > 0) {
    # Divide each column by its corresponding singular value to get orthonormal U
    res$u <- sweep(res$u, 2, res$d, "/")
  }
  
  k_use <- length(res$d)
  if (k_use > 0) {
    pcs <- paste0("PC", seq_len(k_use))
    colnames(res$u) <- pcs
    colnames(res$v) <- pcs
    # Don't set names on d to match multivarious::sdev() behavior
    # names(res$d) <- pcs
  }
  
  res$k <- k_use
  res
}