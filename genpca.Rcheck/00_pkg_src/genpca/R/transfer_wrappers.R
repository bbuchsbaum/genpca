#' Compatibility wrapper for multivarious::transfer
#'
#' Provides a transfer method that accepts `source`/`target` arguments
#' used in the tests. The method forwards to `multivarious`'s
#' transfer method which expects `from`/`to`.
#'
#' @param x A cross_projector object.
#' @param new_data New data to transfer.
#' @param from Source space ("X" or "Y"). Also accepts legacy `source` parameter.
#' @param to Target space ("X" or "Y"). Also accepts legacy `target` parameter.
#' @param source Legacy parameter name for `from`. Deprecated, use `from` instead.
#' @param target Legacy parameter name for `to`. Deprecated, use `to` instead.
#' @param opts Options list passed to multivarious::transfer.
#' @param ... Additional arguments passed through.
#' @return Matrix with transferred data.
#' @importFrom multivarious transfer
#' @export
transfer.cross_projector <- function(x, new_data,
                                     from = NULL,
                                     to = NULL,
                                     source = NULL,
                                     target = NULL,
                                     opts = list(),
                                     ...) {
  
  # Handle legacy parameter names (source/target -> from/to)
  if (!is.null(source) && is.null(from)) {
    from <- source
  }
  if (!is.null(target) && is.null(to)) {
    to <- target
  }
  
  # Set defaults if not provided
  if (is.null(from)) from <- "X"
  if (is.null(to)) to <- "Y"
  
  # Match arguments
  from <- match.arg(from, c("X", "Y"))
  to   <- match.arg(to, c("X", "Y"))
  
  # Call multivarious transfer method directly to avoid method dispatch recursion
  getS3method("transfer", "cross_projector", envir = asNamespace("multivarious"))(
    x, new_data, from = from, to = to, opts = opts, ...
  )
}
