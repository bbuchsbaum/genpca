#' Fast generalized matrix decomposition (dense/sparse dispatch)
#'
#' Computes the generalized SVD of X with row metric Q and column metric R,
#' equivalent to the eigendecomposition used by \code{\link{genpca}} with
#' \code{method = "spectra"}. Uses primal (p <= n) or dual (n < p) formulation
#' to minimize computation, with optional Cholesky caching for repeated calls.
#'
#' @section When is this fast:
#' This implementation is faster than the "eigen" method when:
#' \itemize{
#'   \item \code{k << min(n, p)}: Only top-k eigenvalues needed (uses ARPACK)
#'   \item Repeated calls with same Q or R: Cholesky factors are cached
#'   \item Large matrices where full eigendecomposition is expensive
#' }
#' For small matrices or when \code{k} is close to \code{min(n, p)}, the
#' overhead of iterative methods may make "eigen" faster.
#'
#' @param X numeric matrix (n x p)
#' @param Q,R constraints (weights/metrics) for rows/cols. Must be symmetric
#'   positive (semi-)definite. Can be dense matrices, sparse matrices, or
#'   diagonal matrices.
#' @param k number of components to extract (must be >= 1 and <= min(n, p))
#' @param tol tolerance for filtering near-zero singular values. Default 1e-9.
#' @param maxit maximum iterations (ignored, kept for API compatibility)
#' @param seed random seed (ignored, kept for API compatibility)
#' @param topk logical; use top-k symmetric eigen via ARPACK when available.
#'   Defaults to TRUE. Set to FALSE to force full eigendecomposition.
#' @param cache logical; cache Cholesky factors across calls when constraints
#'   are dense. Defaults to TRUE. Use \code{\link{gmd_clear_cache}} to clear.
#'
#' @return A list with components:
#'   \describe{
#'     \item{u}{n x k matrix of scores (left singular vectors scaled)}
#'     \item{v}{p x k matrix of components (right singular vectors, R-scaled)}
#'     \item{d}{length-k vector of singular values}
#'     \item{k}{number of components returned (may be < requested if rank-deficient)}
#'   }
#'
#' @seealso \code{\link{genpca}} for the high-level interface,
#'   \code{\link{gmd_clear_cache}} to clear Cholesky cache
#' @keywords internal
#' @importFrom methods as is
gmd_fast_cpp <- function(X, Q, R, k, tol = 1e-9, maxit = 1000L, seed = 1234L, topk = TRUE, cache = TRUE) {
  # Input validation
  n <- nrow(X); p <- ncol(X)
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a single positive integer >= 1")
  }

  k <- as.integer(k)
  if (k > min(n, p)) {
    warning("k (", k, ") exceeds min(n, p) = ", min(n, p), "; will return at most ", min(n, p), " components")
    k <- min(n, p)
  }

  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)

  # coerce symmetric dense forms that trip Rcpp
  if (inherits(Q, "dsyMatrix")) Q <- methods::as(Q, "dgeMatrix")
  if (inherits(R, "dsyMatrix")) R <- methods::as(R, "dgeMatrix")

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
  
  # Set k to the number of components found
  res$k <- length(res$d)
  res
}