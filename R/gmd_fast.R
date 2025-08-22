#' Fast generalized matrix decomposition (dense/sparse dispatch wrapper)
#' This overrides the C++ exported function to handle dsyMatrix properly
#' @param X numeric matrix (n x p)
#' @param Q,R constraints (weights/metrics) for rows/cols
#' @param k number of components
#' @param tol tolerance
#' @param maxit maximum iterations (for compatibility, passed to old implementation if needed)
#' @param seed random seed (for compatibility, ignored in new implementation)
#' @keywords internal
#' @importFrom methods as is
gmd_fast_cpp <- function(X, Q, R, k, tol = 1e-9, maxit = 1000L, seed = 1234L) {
  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)
  
  # Coerce symmetric dense forms that trip Rcpp
  if (inherits(Q, "dsyMatrix")) Q <- methods::as(Q, "dgeMatrix")
  if (inherits(R, "dsyMatrix")) R <- methods::as(R, "dgeMatrix")
  
  # Dispatch based on sparsity
  if (methods::is(Q, "sparseMatrix") || methods::is(R, "sparseMatrix")) {
    Q <- methods::as(Q, "dgCMatrix")
    R <- methods::as(R, "dgCMatrix")
    return(gmd_fast_cpp_sp(X, Q, R, k, tol))
  } else {
    return(gmd_fast_cpp_dn(X, as.matrix(Q), as.matrix(R), k, tol))
  }
}