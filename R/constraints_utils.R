#' @title Utilities for constraints
#' @description Helpers to make constraint matrices symmetric positive definite (SPD).
#' @name constraints_utils
#' @keywords internal
NULL

#' @title Check if matrix is SPD
#' @description Check if a matrix is symmetric positive semi-definite (within
#' `tol`) by attempting a Cholesky factorization of `A + tol*scale*I`.
#' @param A numeric matrix or Matrix::Matrix
#' @param tol relative tolerance: eigenvalues above `-tol * max(abs(diag(A)))`
#'   are treated as non-negative, so PSD-but-singular metrics are accepted
#' @return logical TRUE if symmetric PSD within tolerance, FALSE otherwise
#' @keywords internal
is_spd <- function(A, tol = 1e-6) {
  if (!inherits(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = FALSE)
  }
  if (!Matrix::isSymmetric(A)) {
    return(FALSE)
  }
  # Probe A + shift*I: positive definite iff min eigenvalue of A > -shift.
  # The factorization must be one that FAILS on indefinite input. For dense
  # matrices Matrix::Cholesky() dispatches to LAPACK's pivoted dpstrf, which
  # succeeds (with only a warning) on indefinite matrices, so dense input must
  # go through base::chol (dpotrf) instead.
  scale <- max(abs(Matrix::diag(A)), 0)
  if (!is.finite(scale)) return(FALSE)
  if (scale == 0) scale <- 1
  probe <- A + Matrix::Diagonal(nrow(A), x = tol * scale)
  if (methods::is(probe, "sparseMatrix")) {
    # CHOLMOD errors on non-PD input; warnings are near-singular noise.
    suppressWarnings(tryCatch({
      Matrix::Cholesky(probe, LDL = FALSE, Imult = 0, super = TRUE)
      TRUE
    }, error = function(e) FALSE))
  } else {
    tryCatch({
      chol(as.matrix(probe))
      TRUE
    }, error = function(e) FALSE)
  }
}

#' @title Coerce to general CSC sparse matrix (dgCMatrix)
#' @description Replacement for the deprecated direct `as(., "dgCMatrix")`
#' coercion from symmetric/triangular/diagonal Matrix classes.
#' @param A a matrix or Matrix
#' @return a dgCMatrix
#' @keywords internal
as_dgc <- function(A) {
  if (inherits(A, "dgCMatrix")) return(A)
  if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  methods::as(methods::as(methods::as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

#' @title Coerce dense symmetric Matrix classes to general dense (dgeMatrix)
#' @description Replacement for the deprecated direct `as(., "dgeMatrix")`
#' coercion from dsyMatrix/dpoMatrix.
#' @param A a dense Matrix
#' @return a dgeMatrix
#' @keywords internal
as_dge <- function(A) {
  if (inherits(A, "dgeMatrix")) return(A)
  methods::as(methods::as(A, "generalMatrix"), "unpackedMatrix")
}

#' @title Clip a symmetric matrix to the PSD cone
#' @description Spectral clip: eigen-decompose and set negative eigenvalues to
#' zero. Unlike [ensure_spd()] (a diagonal ridge shift), this preserves the
#' non-negative part of the spectrum exactly. Requires a dense
#' eigendecomposition, so large sparse matrices are refused.
#' @param M numeric matrix or Matrix::Matrix
#' @param tol tolerance passed to [is_spd()] for the fast-path check
#' @param dense_maxn refuse sparse input larger than this (clip densifies)
#' @return a dense Matrix, symmetric PSD
#' @keywords internal
clip_psd <- function(M, tol = 1e-6, dense_maxn = 2000L) {
  if (!inherits(M, "Matrix")) {
    M <- Matrix::Matrix(M, sparse = FALSE)
  }
  M <- Matrix::forceSymmetric(M, uplo = "U")
  if (is_spd(M, tol = tol)) return(M)
  n <- nrow(M)
  if (methods::is(M, "sparseMatrix") && n > dense_maxn) {
    stop("constraints_remedy = 'clip' needs a dense eigendecomposition, but ",
         "the matrix is sparse with ", n, " rows (> ", dense_maxn,
         "). Use constraints_remedy = 'ridge' instead.")
  }
  ee <- eigen(as.matrix(M), symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  Mc <- ee$vectors %*% (vals * t(ee$vectors))
  Matrix::Matrix((Mc + t(Mc)) / 2, sparse = FALSE)
}

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
  xvals <- if (methods::is(M, "sparseMatrix")) M@x else as.numeric(M)
  if (length(xvals) && !all(is.finite(xvals))) {
    stop("ensure_spd: matrix contains non-finite values.")
  }
  M <- Matrix::forceSymmetric(M, uplo = "U")

  n <- nrow(M)

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
