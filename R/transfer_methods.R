#' Transfer data between domains using a cross_projector
#'
#' This method wraps [multivarious::transfer()] but allows the argument names
#' `source` and `target` instead of `from` and `to` for convenience.
#'
#' @inherit multivarious::transfer
#' @param object A cross_projector object.
#' @param new_data Matrix or data.frame to transfer.
#' @param source Source domain ("X" or "Y").
#' @param target Target domain ("X" or "Y").
#' @param ... Additional arguments passed through.
#' @export
transfer.cross_projector <- function(object, new_data,
                                     source = c("X", "Y"),
                                     target = c("X", "Y"), ...) {
  source <- match.arg(source)
  target <- match.arg(target)
  multivarious:::transfer.cross_projector(object, new_data,
                                          from = source, to = target, ...)
}
