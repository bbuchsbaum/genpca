#' Generalized PLS via Implicit Operator (PLS-SVD / GPLSSVD)
#'
#' Canonical (two-block) generalized PLS using sparse-friendly implicit
#' matrix-vector products. Solves the SVD of the operator \eqn{S = Xe' Ye}
#' without materializing \eqn{Xe = Mx^{1/2} X Ax^{1/2}} or \eqn{Ye = My^{1/2} Y Ay^{1/2}}.
#'
#' This follows the GPLSSVD/PLS-SVD formulation (Beaton, eqs. 10–14):
#' the top `ncomp` singular triplets of S are computed by iterative SVD
#' on the linear maps v -> S v and u -> S^T u, implemented with metric
#' Cholesky multiplies/solves when possible. Works with dense or sparse
#' `Matrix` inputs and constraint metrics.
#'
#' Returns a `multivarious::cross_projector` with X-/Y-weights (vx, vy)
#' chosen to provide natural projection of new data (`X %*% vx`, `Y %*% vy`).
#' Additional GPLSSVD quantities are attached to the object for access:
#' singular values `d`, generalized weights `p`, `q`, variable scores `fi`, `fj`,
#' and row latent variables `lx`, `ly`.
#'
#' @param X Numeric or Matrix, n x p.
#' @param Y Numeric or Matrix, n x q. Must have same n as `X`.
#' @param Ax Column metric for X (W_X): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param Ay Column metric for Y (W_Y): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param Mx Row metric for X (M_X): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param My Row metric for Y (M_Y): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param ncomp Number of components to extract (rank-k). Default 2.
#' @param preproc_x,preproc_y Optional `multivarious` preprocessors (e.g., `center()`).
#'   Defaults to `multivarious::pass()` (no-op).
#' @param svd_backend Character, one of `"RSpectra"` (default) or `"irlba"` for
#'   iterative SVD. If neither backend is available, a dense fallback is used
#'   for small problems by materializing S.
#' @param svd_opts List of options passed to the SVD backend, e.g., `tol`, `maxitr`.
#' @param verbose Logical; print brief progress messages.
#'
#' @return An object of class `c("genpls", "cross_projector", "projector")` with:
#'   \describe{
#'     \item{vx, vy}{X- and Y- weights usable with predict/transfer (stored in cross_projector)}
#'     \item{d}{singular values (attached field)}
#'     \item{p, q}{generalized weights \eqn{W_X^{-1/2} u}, \eqn{W_Y^{-1/2} v} (attached)}
#'     \item{fi, fj}{variable/component scores \eqn{W_X p D}, \eqn{W_Y q D} (attached)}
#'     \item{lx, ly}{row latent variables \eqn{M_X^{1/2} X W_X p}, \eqn{M_Y^{1/2} Y W_Y q} (attached)}
#'     \item{metrics}{the supplied metrics (attached)}
#'   }
#'
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE) &&
#'     requireNamespace("multivarious", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100; p <- 40; q <- 30
#'   X <- matrix(rnorm(n*p), n, p)
#'   Y <- matrix(rnorm(n*q), n, q)
#'   w <- runif(n); w <- w/sum(w)
#'   Mx <- My <- Matrix::Diagonal(x = w)
#'   fit <- genpls(X, Y, Mx = Mx, My = My, ncomp = 2,
#'                 preproc_x = multivarious::center(),
#'                 preproc_y = multivarious::center())
#'   fit$d  # singular values
#' }
#'
#' @references
#' Beaton, Dougal. Generalized eigen, singular value, and partial least squares
#' decompositions: The GSVD package. (Eqs. 10–14). 2020.
#'
#' @importFrom Matrix Matrix Diagonal crossprod t forceSymmetric Cholesky solve
#' @importFrom RSpectra svds
#' @importFrom multivarious cross_projector fit fit_transform pass
#' @export
genpls <- function(X, Y,
                   Ax = NULL, Ay = NULL,
                   Mx = NULL, My = NULL,
                   ncomp = 2,
                   preproc_x = multivarious::pass(),
                   preproc_y = multivarious::pass(),
                   svd_backend = c("RSpectra", "irlba"),
                   svd_opts = list(tol = 1e-7, maxitr = 1000),
                   verbose = FALSE) {

  svd_backend <- match.arg(svd_backend)
  stopifnot(length(ncomp) == 1L, ncomp >= 1)

  n <- nrow(X)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows.")

  # Preprocess (e.g., centering/scaling) via multivarious, which requires base
  # matrices. Sparse inputs with a no-op preprocessor skip the densification:
  # the pass() prepper is fitted on a placeholder with the right column count
  # and the sparse data flow to the operator untouched.
  is_pass_prepper <- function(p) {
    inherits(p, "prepper") &&
      all(vapply(p$steps, inherits, logical(1), "pass"))
  }
  prep_block <- function(A, preproc) {
    if (inherits(A, "sparseMatrix") && is_pass_prepper(preproc)) {
      ft <- multivarious::fit_transform(preproc, matrix(0, 1L, ncol(A)))
      list(proc = ft$preproc, A = A)
    } else {
      ft <- multivarious::fit_transform(preproc, as.matrix(A))
      list(proc = ft$preproc, A = ft$transformed)
    }
  }
  bx <- prep_block(X, preproc_x)
  by <- prep_block(Y, preproc_y)
  proc_x <- bx$proc
  proc_y <- by$proc
  Xp <- bx$A
  Yp <- by$A

  # Delegate to operator implementation (already memory-safe)
  if (verbose) message("Computing top-", ncomp, " GPLSSVD components via operator (", svd_backend, ") ...")
  op <- gplssvd_op(Xp, Yp,
                   XLW = Mx, YLW = My,
                   XRW = Ax, YRW = Ay,
                   k = ncomp, center = FALSE, scale = FALSE,
                   svd_backend = svd_backend, svd_opts = svd_opts)

  # Derive projection weights for multivarious wrapper: W_X p = Fi D^{-1}
  # (`rep(invd, each = nrow)` scales column j by 1/d[j])
  invd <- ifelse(op$d > 0, 1 / op$d, 0)
  vx <- op$fi * rep(invd, each = nrow(op$fi))
  vy <- op$fj * rep(invd, each = nrow(op$fj))

  obj <- multivarious::cross_projector(
    vx = vx,
    vy = vy,
    preproc_x = proc_x,
    preproc_y = proc_y,
    classes = "genpls"
  )

  # Attach GPLSSVD details for direct access
  obj$d  <- op$d
  obj$p  <- op$p
  obj$q  <- op$q
  obj$fi <- op$fi
  obj$fj <- op$fj
  obj$lx <- op$lx
  obj$ly <- op$ly
  obj$metrics <- list(Ax = Ax, Ay = Ay, Mx = Mx, My = My)
  obj$ncomp   <- ncomp
  obj$backend <- svd_backend

  if (verbose) message("genpls finished.")
  obj
}
