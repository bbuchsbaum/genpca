#' Matrix-free operator for the irlba backend
#'
#' Internal matrix-like object storing dimensions and functions for forward
#' and adjoint products. Its multiplication methods preserve ordinary matrix
#' product orientation for vectors and blocks without storing the matrix.
#'
#' @name genpca_irlba_operator
#' @aliases genpca_irlba_operator-class dim,genpca_irlba_operator-method
#'   %*%,genpca_irlba_operator,ANY-method %*%,ANY,genpca_irlba_operator-method
#' @keywords internal
NULL

# A matrix-free matrix product for irlba's documented S4 interface. Unlike a
# placeholder matrix plus the retired `mult` argument, the object itself owns
# both products. In particular x %*% A returns a row matrix, not A' %*% x.
methods::setClass("genpca_irlba_operator",
                  slots = c(dimensions = "integer", forward = "function",
                            adjoint = "function"))

methods::setMethod("dim", "genpca_irlba_operator",
                   function(x) x@dimensions)

methods::setMethod("%*%", c(x = "genpca_irlba_operator", y = "ANY"),
                   function(x, y) {
                     y <- as.matrix(y)
                     if (nrow(y) != ncol(x)) stop("non-conformable arguments")
                     as.matrix(x@forward(y))
                   })

methods::setMethod("%*%", c(x = "ANY", y = "genpca_irlba_operator"),
                   function(x, y) {
                     x <- if (is.null(dim(x))) matrix(x, nrow = 1L) else as.matrix(x)
                     if (ncol(x) != nrow(y)) stop("non-conformable arguments")
                     t(as.matrix(y@adjoint(t(x))))
                   })

.irlba_operator <- function(opc) {
  methods::new("genpca_irlba_operator", dimensions = as.integer(opc$dims),
               forward = opc$S_mv, adjoint = opc$ST_mv)
}
