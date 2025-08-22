#' Create Weight Operator Function
#'
#' Returns a closure that applies a weight matrix W or its transformations
#' (square root, inverse, or combinations thereof) to vectors/matrices.
#'
#' @param W A weight matrix (SPD) or NULL for identity
#' @param transpose Logical, whether to transpose W before applying
#' @param sqrt Logical, whether to use square root of W
#' @param inverse Logical, whether to use inverse of W
#' @return A function that applies the requested transformation
#' 
#' @keywords internal
#' @importFrom Matrix isDiagonal chol solve Diagonal
as_weight_operator <- function(W, transpose = FALSE, sqrt = FALSE, inverse = FALSE) {
  # Handle NULL case (identity operator)
  if (is.null(W)) {
    return(function(x) x)
  }
  
  # For diagonal matrices, we can optimize
  if (Matrix::isDiagonal(W)) {
    d <- diag(W)
    
    if (sqrt && inverse) {
      # W^{-1/2}
      d_op <- 1 / sqrt(d)
      return(function(x) {
        if (is.matrix(x) || inherits(x, "Matrix")) {
          # Column-wise multiplication for matrices
          sweep(x, 1, d_op, "*")
        } else {
          x * d_op
        }
      })
    } else if (sqrt && !inverse) {
      # W^{1/2}
      d_op <- sqrt(d)
      return(function(x) {
        if (is.matrix(x) || inherits(x, "Matrix")) {
          sweep(x, 1, d_op, "*")
        } else {
          x * d_op
        }
      })
    } else if (!sqrt && inverse) {
      # W^{-1}
      d_op <- 1 / d
      return(function(x) {
        if (is.matrix(x) || inherits(x, "Matrix")) {
          sweep(x, 1, d_op, "*")
        } else {
          x * d_op
        }
      })
    } else {
      # W
      return(function(x) {
        if (is.matrix(x) || inherits(x, "Matrix")) {
          sweep(x, 1, d, "*")
        } else {
          x * d
        }
      })
    }
  }
  
  # For general matrices
  if (sqrt && inverse) {
    # W^{-1/2} via Cholesky decomposition
    # W = R^T R (where R is upper triangular from chol), so W^{-1/2} = R^{-1}
    W_chol <- Matrix::chol(W)  # Returns upper triangular R where W = R^T R
    W_op <- Matrix::solve(W_chol)  # R^{-1}, so W^{-1/2} = R^{-1}
    if (transpose) W_op <- Matrix::t(W_op)
  } else if (sqrt && !inverse) {
    # W^{1/2} via Cholesky decomposition
    # W = R^T R, so W^{1/2} = R^T
    W_chol <- Matrix::chol(W)
    W_op <- Matrix::t(W_chol)
    if (transpose) W_op <- Matrix::t(W_op)
  } else if (!sqrt && inverse) {
    # W^{-1}
    W_op <- Matrix::solve(W)
    if (transpose) W_op <- Matrix::t(W_op)
  } else {
    # W
    W_op <- W
    if (transpose) W_op <- Matrix::t(W_op)
  }
  
  # Return the operator function
  function(x) {
    result <- W_op %*% x
    # Ensure numeric output for RSpectra
    if (is.vector(x)) {
      as.numeric(result)
    } else {
      result
    }
  }
}