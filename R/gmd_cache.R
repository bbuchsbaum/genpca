#' @title Internal cache for GMD factorizations
#' @keywords internal
.gmd_cache <- new.env(parent = emptyenv())

#' @keywords internal
.digest_dense_matrix <- function(M) {
  # Round to stabilize digest against tiny numeric jitter
  digest::digest(list(dim = dim(M), data = round(as.numeric(M), 8)))
}

#' @keywords internal
.digest_sparse_matrix <- function(M) {
  # Use slots to avoid densifying
  stopifnot(inherits(M, "sparseMatrix"))
  digest::digest(list(dim = M@Dim, i = M@i, p = M@p, x = round(M@x, 8)))
}

#' @title Get (and cache) a *lower* Cholesky factor for a dense SPD matrix
#' @param A numeric or dense Matrix (SPD). If sparse, falls back to dense.
#' @return a base numeric matrix L (lower triangular) with A = L %*% t(L)
#' @keywords internal
get_chol_lower_dense <- function(A) {
  if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = FALSE)
  if (methods::is(A, "sparseMatrix")) A <- methods::as(A, "denseMatrix")
  key <- paste0("L_", .digest_dense_matrix(A))
  if (exists(key, envir = .gmd_cache, inherits = FALSE)) {
    return(get(key, envir = .gmd_cache, inherits = FALSE))
  }
  L <- t(chol(as.matrix(A)))  # base::chol returns upper by default
  assign(key, L, envir = .gmd_cache)
  L
}

#' Clear internal cache for matrix decompositions
#'
#' @description Clears the internal cache used by generalized matrix decomposition functions.
#' This can be useful to free up memory or when working with different datasets.
#'
#' @return Invisibly returns TRUE after clearing the cache.
#'
#' @examples
#' # Clear the internal cache
#' gmd_clear_cache()
#'
#' @export
gmd_clear_cache <- function() {
  rm(list = ls(envir = .gmd_cache, all.names = TRUE), envir = .gmd_cache)
  invisible(TRUE)
}