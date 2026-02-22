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
#' @param maxit maximum iterations for the Spectra eigensolver. Default 1000.
#' @param seed random seed (ignored, kept for API compatibility)
#' @param topk logical; use top-k symmetric eigen via ARPACK when available.
#'   Defaults to TRUE. Set to FALSE to force full eigendecomposition.
#' @param cache logical; cache Cholesky factors across calls when constraints
#'   are dense. Defaults to TRUE. Use \code{\link{gmd_clear_cache}} to clear.
#' @param auto_topk logical; when TRUE (default), use top-k only when
#'   \code{k/min(n,p)} is small and \code{min(n,p)} is large enough.
#' @param topk_ratio threshold used by \code{auto_topk}. If
#'   \code{k/min(n,p) <= topk_ratio}, top-k is used. Default 0.08.
#' @param topk_min_dim minimum \code{min(n,p)} required before top-k is used
#'   under \code{auto_topk}. Default 200.
#' @param diag_fast logical; if TRUE (default) and both constraints are
#'   diagonal, use a weighted-SVD fast path.
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
is_diagonal_metric <- function(M) {
  inherits(M, "ddiMatrix") || Matrix::isDiagonal(M)
}

diag_metric_values <- function(M) {
  as.numeric(Matrix::diag(M))
}

should_use_topk <- function(k, min_dim, topk, auto_topk, topk_ratio, topk_min_dim) {
  if (!isTRUE(topk)) return(FALSE)
  if (!isTRUE(auto_topk)) return(TRUE)
  if (min_dim < topk_min_dim) return(FALSE)
  (k / min_dim) <= topk_ratio
}

gmd_fast_diag <- function(X, q_diag, r_diag, k, tol, maxit, topk) {
  n <- nrow(X)
  p <- ncol(X)
  min_dim <- min(n, p)
  k_use <- min(k, min_dim)

  q_sqrt <- sqrt(pmax(q_diag, 0))
  r_sqrt <- sqrt(pmax(r_diag, 0))
  r_invsqrt <- ifelse(r_sqrt > tol, 1 / r_sqrt, 0)

  # Mirror the C++ primal/dual formulation exactly.
  Xw <- sweep(X, 1, q_sqrt, `*`)
  if (p <= n) {
    Xw <- sweep(Xw, 2, r_invsqrt, `*`)  # primal: Q^{1/2} X R^{-1/2}
  } else {
    Xw <- sweep(Xw, 2, r_sqrt, `*`)     # dual:   Q^{1/2} X R^{1/2}
  }

  sv <- NULL
  if (isTRUE(topk) && k_use < min_dim) {
    sv <- tryCatch(
      RSpectra::svds(Xw, k = k_use, opts = list(maxitr = maxit, tol = tol)),
      error = function(e) NULL
    )
  }

  if (is.null(sv)) {
    sv_full <- base::svd(Xw, nu = k_use, nv = k_use)
    d <- sv_full$d[seq_len(k_use)]
    Vw <- sv_full$v[, seq_len(k_use), drop = FALSE]
  } else {
    d <- sv$d
    Vw <- sv$v
  }

  keep <- d > tol
  if (!any(keep)) {
    return(list(
      u = matrix(0, nrow = n, ncol = 0),
      v = matrix(0, nrow = p, ncol = 0),
      d = numeric(0)
    ))
  }

  d <- as.numeric(d[keep])
  Vw <- Vw[, keep, drop = FALSE]

  # Components in original space: C = R^{1/2} * right-singular-vectors
  components <- sweep(Vw, 1, r_sqrt, `*`)
  # Scores in returned parametrization: U = Q X C
  scores <- sweep(X %*% components, 1, q_diag, `*`)

  list(u = scores, v = components, d = d)
}

gmd_fast_cpp <- function(X, Q, R, k, tol = 1e-9, maxit = 1000L, seed = 1234L,
                         topk = TRUE, cache = TRUE, auto_topk = TRUE,
                         topk_ratio = 0.08, topk_min_dim = 200L,
                         diag_fast = TRUE) {
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
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1 || maxit != floor(maxit)) {
    stop("maxit must be a single positive integer >= 1")
  }
  maxit <- as.integer(maxit)
  if (!is.numeric(topk_ratio) || length(topk_ratio) != 1 || topk_ratio <= 0 || topk_ratio > 1) {
    stop("topk_ratio must be a single number in (0, 1].")
  }
  if (!is.numeric(topk_min_dim) || length(topk_min_dim) != 1 || topk_min_dim < 2 || topk_min_dim != floor(topk_min_dim)) {
    stop("topk_min_dim must be a single integer >= 2.")
  }
  topk_min_dim <- as.integer(topk_min_dim)

  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)

  # coerce symmetric dense forms that trip Rcpp
  if (inherits(Q, "dsyMatrix")) Q <- methods::as(Q, "dgeMatrix")
  if (inherits(R, "dsyMatrix")) R <- methods::as(R, "dgeMatrix")

  min_dim <- min(n, p)
  use_topk <- should_use_topk(k, min_dim, topk, auto_topk, topk_ratio, topk_min_dim)

  if (isTRUE(diag_fast) && is_diagonal_metric(Q) && is_diagonal_metric(R)) {
    res <- gmd_fast_diag(
      X = X,
      q_diag = diag_metric_values(Q),
      r_diag = diag_metric_values(R),
      k = k,
      tol = tol,
      maxit = maxit,
      topk = use_topk
    )
  } else {
  primal <- (p <= n)  # use primal when small side is p

  # ---- choose path and apply caching for dense constraints ----
  if (primal) {
    if (!methods::is(R, "sparseMatrix") && isTRUE(cache)) {
      L_R <- get_chol_lower_dense(R)
      if (methods::is(Q, "sparseMatrix")) {
        res <- gmd_fast_cpp_primal_sp(X, methods::as(Q, "dgCMatrix"), L_R, k, tol, maxit, use_topk)
      } else {
        res <- gmd_fast_cpp_primal_dn(X, as.matrix(Q), L_R, k, tol, maxit, use_topk)
      }
    } else {
      # fall back to non-cached path
      if (methods::is(Q, "sparseMatrix") || methods::is(R, "sparseMatrix")) {
        res <- gmd_fast_cpp_sp(X, methods::as(Q, "dgCMatrix"), methods::as(R, "dgCMatrix"), k, tol, maxit, use_topk)
      } else {
        res <- gmd_fast_cpp_dn(X, as.matrix(Q), as.matrix(R), k, tol, maxit, use_topk)
      }
    }
  } else {
    # dual path (n < p): cache Cholesky for Q when dense
    if (!methods::is(Q, "sparseMatrix") && isTRUE(cache)) {
      L_Q <- get_chol_lower_dense(Q)
      if (methods::is(R, "sparseMatrix")) {
        res <- gmd_fast_cpp_dual_sp(X, L_Q, methods::as(R, "dgCMatrix"), k, tol, maxit, use_topk)
      } else {
        res <- gmd_fast_cpp_dual_dn(X, L_Q, as.matrix(R), k, tol, maxit, use_topk)
      }
    } else {
      if (methods::is(Q, "sparseMatrix") || methods::is(R, "sparseMatrix")) {
        res <- gmd_fast_cpp_sp(X, methods::as(Q, "dgCMatrix"), methods::as(R, "dgCMatrix"), k, tol, maxit, use_topk)
      } else {
        res <- gmd_fast_cpp_dn(X, as.matrix(Q), as.matrix(R), k, tol, maxit, use_topk)
      }
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
