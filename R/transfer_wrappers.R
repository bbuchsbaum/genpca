#' Compatibility wrapper for multivarious::transfer
#'
#' Provides a transfer method that accepts `source`/`target` arguments
#' used in the tests. The method forwards to `multivarious`'s
#' `transfer.cross_projector()` which expects `from`/`to`.
#'
#' @param object A cross_projector object.
#' @param newdata New data to transfer.
#' @param source Source space ("X" or "Y").
#' @param target Target space ("X" or "Y").
#' @param ... Additional arguments passed through.
#' @return Matrix with transferred data.
#' @export
transfer.cross_projector <- function(object, newdata,
                                     source = c("X", "Y"),
                                     target = c("Y", "X"),
                                     ...) {
  from <- match.arg(source)
  to   <- match.arg(target)
  multivarious:::transfer.cross_projector(object = object,
                                          newdata = newdata,
                                          from = from,
                                          to = to,
                                          ...)
}
