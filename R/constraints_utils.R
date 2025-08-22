#' @title Utilities for constraints
#' @description Helpers to make constraint matrices symmetric positive definite (SPD).
#' @keywords internal
NULL

#' @title Ensure SPD (sparse-friendly)
#' @description Force a symmetric matrix to be symmetric positive definite (SPD).
#' Uses a Gershgorin-based diagonal shift; falls back to nearPD for small dense matrices.
#' @param M numeric matrix or Matrix::Matrix
#' @param tol jitter tolerance (default 1e-8)
#' @param nearpd_maxn only use nearPD when n <= nearpd_maxn and matrix is dense
#' @return a Matrix object (sparse stays sparse when possible)
#' @keywords internal
#' @importFrom Matrix forceSymmetric Diagonal rowSums Cholesky nearPD
#' @importFrom methods as is
ensure_spd <- function(M, tol = 1e-6, nearpd_maxn = 2000L) {
  if (!inherits(M, "Matrix")) {
    M <- Matrix::Matrix(M, sparse = FALSE)
  }
  M <- Matrix::forceSymmetric(M, uplo = "U")
  
  n <- nrow(M)
  
  # Quick SPD probe via Cholesky (fast if already SPD)
  is_spd <- function(A) {
    ok <- TRUE
    tryCatch({
      Matrix::Cholesky(A, LDL = FALSE, Imult = 0, super = TRUE)
    }, error = function(e) ok <<- FALSE)
    ok
  }
  
  if (is_spd(M)) return(M)
  
  # Gershgorin-based diagonal shift to guarantee SPD without eigen computations
  # Gershgorin circle theorem: eigenvalues lie in union of circles centered at a_ii 
  # with radius sum_j≠i |a_ij|
  # Making a_ii - sum_j≠i |a_ij| > 0 ensures strict diagonal dominance => SPD
  d <- Matrix::diag(M)
  rs <- Matrix::rowSums(abs(M)) - abs(d)
  min_margin <- suppressWarnings(min(d - rs))
  
  if (!is.finite(min_margin)) min_margin <- -1
  
  if (min_margin <= 0) {
    # Use a more conservative shift for better numerical stability
    shift <- (-min_margin) + max(tol, abs(min_margin) * 0.1)
    M <- M + Matrix::Diagonal(n, x = shift)
  }
  
  if (is_spd(M)) return(M)
  
  # Final fallback for small dense matrices
  if (!methods::is(M, "sparseMatrix") && n <= nearpd_maxn) {
    Mnp <- Matrix::nearPD(as.matrix(M), corr = FALSE, keepDiag = TRUE)$mat
    return(Matrix::Matrix(Mnp, sparse = FALSE))
  }
  
  # Escalating jitter with more aggressive increases
  jitter <- max(tol, 1e-10 * max(abs(Matrix::diag(M))))
  for (i in seq_len(10)) {
    M2 <- M + Matrix::Diagonal(n, x = jitter)
    if (is_spd(M2)) return(M2)
    jitter <- jitter * 100  # More aggressive scaling
  }
  
  stop("ensure_spd: Unable to make matrix SPD; consider revising constraint.")
}