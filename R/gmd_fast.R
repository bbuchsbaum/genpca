#' Check whether a metric matrix is diagonal
#' @param M a matrix
#' @return logical
#' @keywords internal
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
  q_invsqrt <- ifelse(q_sqrt > tol, 1 / q_sqrt, 0)
  r_invsqrt <- ifelse(r_sqrt > tol, 1 / r_sqrt, 0)

  # Match gmdLA target: singular values of Q^{1/2} X R^{1/2}.
  Xw <- sweep(X, 1, q_sqrt, `*`)
  Xw <- sweep(Xw, 2, r_sqrt, `*`)

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
    Uw <- sv_full$u[, seq_len(k_use), drop = FALSE]
    Vw <- sv_full$v[, seq_len(k_use), drop = FALSE]
  } else {
    d <- sv$d
    Uw <- sv$u
    Vw <- sv$v
  }

  keep <- d > tol
  if (!any(keep)) {
    return(list(
      u = matrix(0, nrow = n, ncol = 0),
      v = matrix(0, nrow = p, ncol = 0),
      ou = matrix(0, nrow = n, ncol = 0),
      ov = matrix(0, nrow = p, ncol = 0),
      d = numeric(0)
    ))
  }

  d <- as.numeric(d[keep])
  Uw <- Uw[, keep, drop = FALSE]
  Vw <- Vw[, keep, drop = FALSE]

  # Metric-orthonormal factors.
  ov <- sweep(Vw, 1, r_invsqrt, `*`)
  ou <- sweep(Uw, 1, q_invsqrt, `*`)
  # Components in original space: C = R^{1/2} * right-singular-vectors
  components <- sweep(Vw, 1, r_sqrt, `*`)
  # Scores in returned parametrization: U = Q X C
  scores <- sweep(X %*% components, 1, q_diag, `*`)

  list(u = scores, v = components, ou = ou, ov = ov, d = d)
}

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
gmd_fast_cpp <- function(X, Q, R, k, tol = 1e-9, maxit = 1000L, seed = 1234L,
                         topk = TRUE, cache = TRUE, auto_topk = TRUE,
                         topk_ratio = 0.08, topk_min_dim = 200L,
                         diag_fast = TRUE) {
  # Input validation
  n <- nrow(X)
  p <- ncol(X)
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
    if (!is.null(res$ou)) res$ou <- res$ou[, keep, drop = FALSE]
    if (!is.null(res$ov)) res$ov <- res$ov[, keep, drop = FALSE]
    res$d <- res$d[keep]
  }

  # Set k to the number of components found
  res$k <- length(res$d)
  res
}

# Metric orthonormalization for small block matrices.
metric_orthonormalize <- function(A, applyM, jitter = 1e-10, tol = 1e-12) {
  if (ncol(A) == 0) return(A)

  MA <- applyM(A)
  G <- Matrix::crossprod(A, MA)
  G <- 0.5 * (G + t(G))
  G <- as.matrix(G)

  # Fast path: Cholesky on the small Gram matrix.
  jitter_now <- max(jitter, 0)
  for (i in 0:6) {
    G_reg <- G
    if (jitter_now > 0) {
      diag(G_reg) <- diag(G_reg) + jitter_now
    }
    C <- try(chol(G_reg), silent = TRUE)
    if (!inherits(C, "try-error")) {
      return(A %*% solve(C))
    }
    jitter_now <- if (jitter_now > 0) jitter_now * 10 else 1e-10
  }

  # Robust fallback for semidefinite / rank-deficient blocks.
  eg <- eigen(G, symmetric = TRUE)
  keep <- which(eg$values > tol)
  if (length(keep) == 0) {
    return(matrix(0.0, nrow = nrow(A), ncol = 0))
  }
  B <- eg$vectors[, keep, drop = FALSE]
  B <- sweep(B, 2, sqrt(eg$values[keep]), "/")
  A %*% B
}

random_sign_matrix <- function(nrow, ncol, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(sample(c(-1, 1), size = nrow * ncol, replace = TRUE), nrow = nrow, ncol = ncol)
}

gmd_randomized_polish <- function(X, applyQ, applyR, U, V, iters = 1L, jitter = 1e-10,
                                  polish_tol = 0) {
  if (iters <= 0L || ncol(U) == 0) {
    return(list(u = U, v = V, d = numeric(ncol(U))))
  }

  d <- numeric(ncol(U))
  d_prev <- NULL
  for (it in seq_len(iters)) {
    Y <- X %*% applyR(V)
    U <- metric_orthonormalize(Y, applyQ, jitter = jitter)
    Z <- Matrix::crossprod(X, applyQ(U))
    V <- metric_orthonormalize(Z, applyR, jitter = jitter)
    RV <- applyR(V)
    # T = U^T Q X R V; reuse Z = X^T Q U to avoid another pass through X.
    T <- Matrix::crossprod(Z, RV)
    sv <- svd(as.matrix(T))
    U <- U %*% sv$u
    V <- V %*% sv$v
    d <- as.numeric(sv$d)

    if (!is.null(d_prev) && polish_tol > 0 && length(d_prev) == length(d)) {
      rel_change <- max(abs(d - d_prev) / pmax(abs(d_prev), 1e-12))
      if (is.finite(rel_change) && rel_change < polish_tol) {
        break
      }
    }
    d_prev <- d
  }

  list(u = U, v = V, d = d)
}

# Randomized low-pass GMD backend in pure R (fallback):
# 2 passes when n_power = 0, 2 + 2*n_power passes otherwise.
gmd_randomized_r <- function(X, Q, R, k,
                             oversample = 20L,
                             n_power = 1L,
                             n_polish = 0L,
                             jitter = 1e-10,
                             tol = 1e-9,
                             polish_tol = 0,
                             seed = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
    stop("k must be a positive integer.")
  }
  if (!is.numeric(oversample) || length(oversample) != 1 || oversample < 0 || oversample != floor(oversample)) {
    stop("oversample must be a non-negative integer.")
  }
  if (!is.numeric(n_power) || length(n_power) != 1 || n_power < 0 || n_power != floor(n_power)) {
    stop("n_power must be a non-negative integer.")
  }
  if (!is.numeric(n_polish) || length(n_polish) != 1 || n_polish < 0 || n_polish != floor(n_polish)) {
    stop("n_polish must be a non-negative integer.")
  }
  if (!is.numeric(polish_tol) || length(polish_tol) != 1 || polish_tol < 0) {
    stop("polish_tol must be a single non-negative number.")
  }

  k <- as.integer(min(k, n, p))
  ell <- as.integer(min(n, p, k + oversample))
  if (ell < 1) {
    return(list(
      u = matrix(0.0, nrow = n, ncol = 0),
      v = matrix(0.0, nrow = p, ncol = 0),
      d = numeric(0),
      k = 0L
    ))
  }

  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)
  applyQ <- function(B) Q %*% B
  applyR <- function(B) R %*% B

  Omega <- random_sign_matrix(p, ell, seed = seed)
  Y <- X %*% applyR(Omega)

  if (n_power > 0L) {
    for (it in seq_len(as.integer(n_power))) {
      Utmp <- metric_orthonormalize(Y, applyQ, jitter = jitter)
      Z <- Matrix::crossprod(X, applyQ(Utmp))
      Y <- X %*% applyR(Z)
    }
  }

  U0 <- metric_orthonormalize(Y, applyQ, jitter = jitter)
  if (ncol(U0) == 0) {
    return(list(
      u = matrix(0.0, nrow = n, ncol = 0),
      v = matrix(0.0, nrow = p, ncol = 0),
      d = numeric(0),
      k = 0L
    ))
  }

  B <- Matrix::crossprod(X, applyQ(U0))
  RB <- applyR(B)
  G <- Matrix::crossprod(B, RB)
  G <- 0.5 * (G + t(G))

  eg <- eigen(as.matrix(G), symmetric = TRUE)
  k_use <- min(k, ncol(eg$vectors))
  if (k_use < 1) {
    return(list(
      u = matrix(0.0, nrow = n, ncol = 0),
      v = matrix(0.0, nrow = p, ncol = 0),
      d = numeric(0),
      k = 0L
    ))
  }

  S <- eg$vectors[, seq_len(k_use), drop = FALSE]
  d <- sqrt(pmax(eg$values[seq_len(k_use)], 0))

  U <- U0 %*% S
  V <- B %*% S
  nz <- d > tol
  if (any(nz)) {
    V[, nz] <- sweep(V[, nz, drop = FALSE], 2, d[nz], "/")
  }
  if (any(!nz)) {
    V[, !nz] <- 0.0
  }

  if (n_polish > 0L) {
    pol <- gmd_randomized_polish(
      X = X,
      applyQ = applyQ,
      applyR = applyR,
      U = U,
      V = V,
      iters = as.integer(n_polish),
      jitter = jitter,
      polish_tol = polish_tol
    )
    U <- pol$u
    V <- pol$v
    d <- pol$d
  }

  keep <- d > tol
  if (!any(keep)) {
    return(list(
      u = matrix(0.0, nrow = n, ncol = 0),
      v = matrix(0.0, nrow = p, ncol = 0),
      d = numeric(0),
      k = 0L
    ))
  }

  U <- U[, keep, drop = FALSE]
  V <- V[, keep, drop = FALSE]
  d <- as.numeric(d[keep])

  list(u = U, v = V, d = d, k = length(d))
}

gmd_randomized <- function(X, Q, R, k,
                           oversample = 20L,
                           n_power = 1L,
                           n_polish = 0L,
                           jitter = 1e-10,
                           tol = 1e-9,
                           polish_tol = 0,
                           seed = NULL,
                           use_cpp = TRUE) {
  if (isTRUE(use_cpp) && exists("gmd_randomized_cpp_dn", mode = "function")) {
    if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
    if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)

    seed_val <- if (is.null(seed)) 1234L else as.integer(seed)
    cpp_res <- tryCatch({
      if (methods::is(Q, "sparseMatrix") && methods::is(R, "sparseMatrix")) {
        gmd_randomized_cpp_sp(
          X = X,
          Q = methods::as(Q, "dgCMatrix"),
          R = methods::as(R, "dgCMatrix"),
          k = as.integer(k),
          oversample = as.integer(oversample),
          n_power = as.integer(n_power),
          n_polish = as.integer(n_polish),
          jitter = jitter,
          tol = tol,
          polish_tol = polish_tol,
          seed = seed_val
        )
      } else if (methods::is(Q, "sparseMatrix")) {
        gmd_randomized_cpp_qsp_rdn(
          X = X,
          Q = methods::as(Q, "dgCMatrix"),
          R = as.matrix(R),
          k = as.integer(k),
          oversample = as.integer(oversample),
          n_power = as.integer(n_power),
          n_polish = as.integer(n_polish),
          jitter = jitter,
          tol = tol,
          polish_tol = polish_tol,
          seed = seed_val
        )
      } else if (methods::is(R, "sparseMatrix")) {
        gmd_randomized_cpp_qdn_rsp(
          X = X,
          Q = as.matrix(Q),
          R = methods::as(R, "dgCMatrix"),
          k = as.integer(k),
          oversample = as.integer(oversample),
          n_power = as.integer(n_power),
          n_polish = as.integer(n_polish),
          jitter = jitter,
          tol = tol,
          polish_tol = polish_tol,
          seed = seed_val
        )
      } else {
        gmd_randomized_cpp_dn(
          X = X,
          Q = as.matrix(Q),
          R = as.matrix(R),
          k = as.integer(k),
          oversample = as.integer(oversample),
          n_power = as.integer(n_power),
          n_polish = as.integer(n_polish),
          jitter = jitter,
          tol = tol,
          polish_tol = polish_tol,
          seed = seed_val
        )
      }
    }, error = function(e) {
      warning("C++ randomized backend failed, falling back to R implementation: ", e$message)
      NULL
    })

    if (!is.null(cpp_res)) {
      if (is.matrix(cpp_res$d)) cpp_res$d <- as.vector(cpp_res$d)
      cpp_res$k <- length(cpp_res$d)
      return(cpp_res)
    }
  }

  gmd_randomized_r(
    X = X,
    Q = Q,
    R = R,
    k = k,
    oversample = oversample,
    n_power = n_power,
    n_polish = n_polish,
    jitter = jitter,
    tol = tol,
    polish_tol = polish_tol,
    seed = seed
  )
}
