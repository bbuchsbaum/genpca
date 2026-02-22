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
  
  # For general symmetric PSD matrices use shared metric operators
  ops <- .metric_operators(W)
  opfun <- if (sqrt && inverse) {
    ops$mult_invsqrt
  } else if (sqrt && !inverse) {
    ops$mult_sqrt
  } else if (!sqrt && inverse) {
    function(x) Matrix::solve(W, x)
  } else {
    ops$mult
  }

  # Transpose is a no-op for symmetric PSD; include for API compatibility
  function(x) opfun(x)
}
