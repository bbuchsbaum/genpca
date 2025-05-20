#' Generalized Partial Least Squares "Correlation" (All Factors, Closed Form) with Adaptive Rank
#'
#' This function computes a \emph{Generalized PLS Correlation} decomposition
#' of two data blocks, \eqn{X} and \eqn{Y}, each with optional row/column
#' constraints: \eqn{M_x, A_x} for \code{X}, and \eqn{M_y, A_y} for \code{Y}.
#' The approach is analogous to the classical \emph{PLS correlation} method,
#' but in an embedding induced by the constraints. Specifically, we define:
#' \deqn{
#'   X^* = \sqrt{M_x}\,X\,\sqrt{A_x},
#'   \quad
#'   Y^* = \sqrt{M_y}\,Y\,\sqrt{A_y},
#' }
#' then compute the cross-product \eqn{\mathbf{C} = (X^*)^\top \,Y^*}, perform an
#' SVD on \eqn{\mathbf{C}}, and map results back to the original domain.
#'
#' In contrast to iterative methods, this yields all factors simultaneously
#' (up to the rank of \(\mathbf{C}\)). For large \eqn{p} or \eqn{q}, building
#' \(\mathbf{C}\) of dimension \eqn{p \times q} can be memory-intensive. One
#' may consider a partial SVD or operator-based approach to avoid forming
#' \(\mathbf{C}\) explicitly. This function demonstrates the direct
#' closed-form approach for moderate sizes.
#'
#' @section Adaptive Rank:
#' For each constraint matrix (\code{Mx}, \code{Ax}, \code{My}, \code{Ay}), if its corresponding
#' \code{rank_*} parameter is \emph{not} a positive integer (i.e., it is \code{NULL}, \code{NA},
#' or \code{0}), an \emph{adaptive} partial eigen approach is used. This approach computes
#' eigenpairs up to \code{max_k} and then keeps only as many as needed to capture
#' \code{var_threshold} fraction of the total variance. If \code{rank_*} is a positive integer,
#' that exact rank is used.
#'
#' @param X A numeric matrix of shape \eqn{(n \times p)}. Must be numeric (not a data.frame).
#' @param Y A numeric matrix of shape \eqn{(n \times q)}. Must have the same \code{n} as \code{X}.
#'
#' @param Mx An \eqn{(n \times n)} PSD matrix for row constraints on \code{X}.
#'        If \code{NULL}, defaults to identity. Similarly for \code{Ax}, \code{My}, and \code{Ay}.
#' @param Ax An \eqn{(p \times p)} PSD matrix for column constraints on \code{X}.
#' @param My An \eqn{(n \times n)} PSD matrix for row constraints on \code{Y}.
#' @param Ay An \eqn{(q \times q)} PSD matrix for column constraints on \code{Y}.
#'
#' @param rank_Mx, rank_Ax, rank_My, rank_Ay 
#'   Either a positive integer (explicit rank) or \code{NULL}/\code{NA}/\code{0} for adaptive approach.
#'   Defaults to \code{NULL}, meaning we do an adaptive method capturing \code{var_threshold}
#'   fraction of variance, capped by \code{max_k}.
#'
#' @param var_threshold Numeric fraction of total variance to capture if ranks are \code{NULL}/\code{NA}/\code{0}.
#'        Defaults to 0.99 (99\%).
#' @param max_k Integer maximum cap for the partial eigen expansions. Defaults to 200.
#'
#' @param ncomp Number of factors (SVD components) to return, up to
#'        \(\min(p,q)\). If \code{ncomp <= 0}, returns all components.
#'
#' @param preproc_x, preproc_y Pre-processing for \code{X} and \code{Y},
#'        e.g. \code{\link[multivarious]{center}}, or \code{\link[multivarious]{pass}}
#'        for no transform. Typically columns are centered for PLS.
#'
#' @param tol Tolerance for small singular values; used to ignore near-zero
#'        eigen/singular values. Defaults to \code{1e-9}.
#' @param maxiter An integer. Ignored in this function but included for interface
#'        consistency with \code{genpls()}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @details 
#' We try to detect if a constraint matrix is identity or diagonal. If it is,
#' we skip partial eigen expansions. Otherwise, we approximate \eqn{\sqrt{M}}
#' via partial eigen. The rank is chosen adaptively or from a user-supplied integer.
#'
#' Then we row-transform \(\mathbf{X}\) with \(\sqrt{M_x}\) and column-transform
#' with \(\sqrt{A_x}\), and do similarly for \(\mathbf{Y}\). We form the cross-product
#' \(\mathbf{C} = (X^*)^\top Y^*\), do an SVD (up to \code{ncomp}), and map loadings
#' back to the original domain.
#'
#' @return A \code{cross_projector} object of class
#'   \code{c("genplscor","cross_projector","projector")} with typical fields:
#'   \describe{
#'     \item{\code{vx}}{A \eqn{(p \times ncomp)} matrix of X-loadings in original space.}
#'     \item{\code{vy}}{A \eqn{(q \times ncomp)} matrix of Y-loadings in original space.}
#'     \item{\code{tilde_vx, tilde_vy}}{The loadings in the embedded domain (\(U, V\)).}
#'     \item{\code{Tx, Ty}}{Optionally, the factor scores in the embedded domain.}
#'     \item{\code{ncomp, rank_Mx, ...}}{The final number of components plus rank details.}
#'     \item{\code{preproc_x, preproc_y}}{The \pkg{multivarious} pre-processor objects.}
#'   }
#'
#' @section Memory & Large Data:
#' This function explicitly builds \(\mathbf{C}\) of size \eqn{p \times q}. For large
#' \(\,(p,q)\), that can be huge. A partial SVD or operator-based approach would be more
#' memory-efficient, but is not shown here.
#'
#' @seealso \code{\link{genpls}} for the iterative approach with deflation + Gram–Schmidt orthonormalization.
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   n <- 10; p <- 5; q <- 4
#'   X <- matrix(rnorm(n*p), n, p)
#'   Y <- matrix(rnorm(n*q), n, q)
#'
#'   # rank_* = NULL => adaptive partial eigen
#'   fit_corr <- genplscorr(X, Y, ncomp=2, verbose=TRUE)
#'   str(fit_corr)
#' }
#'
#' @importFrom Matrix crossprod Diagonal
#' @importFrom multivarious prep pass init_transform
#' @export
genplscorr <- function(X, Y,
                       Mx=NULL, Ax=NULL,
                       My=NULL, Ay=NULL,
                       ncomp=2,
                       preproc_x=pass(),
                       preproc_y=pass(),
                       rank_Mx=NULL, rank_Ax=NULL,
                       rank_My=NULL, rank_Ay=NULL,
                       var_threshold=0.99,
                       max_k=200,
                       tol=1e-9,
                       maxiter=200,  # not actually used
                       verbose=FALSE)
{
  ### 1) Basic checks
  if(!is.matrix(X) || !is.numeric(X)) {
    stop("genplscorr: X must be a numeric matrix.")
  }
  if(!is.matrix(Y) || !is.numeric(Y)) {
    stop("genplscorr: Y must be a numeric matrix.")
  }
  if(nrow(X) != nrow(Y)) {
    stop("genplscorr: X,Y must have the same number of rows.")
  }
  if(ncomp <= 0) {
    ncomp <- min(ncol(X), ncol(Y))
    if(verbose) message("genplscorr: ncomp <= 0 => using all comps up to min(p,q).")
  }
  
  ### 2) Preprocessing
  px <- prep(preproc_x)
  py <- prep(preproc_y)
  Xp <- init_transform(px, X)
  Yp <- init_transform(py, Y)
  
  n <- nrow(Xp)
  p <- ncol(Xp)
  q <- ncol(Yp)
  
  ### 3) Default constraints => identity
  if(is.null(Mx)) Mx <- Matrix::Diagonal(n)
  if(is.null(Ax)) Ax <- Matrix::Diagonal(p)
  if(is.null(My)) My <- Matrix::Diagonal(n)
  if(is.null(Ay)) Ay <- Matrix::Diagonal(q)
  
  ### 4) Use the shared utilities from plsutils.R
  #    build_sqrt_mult, build_invsqrt_mult, etc.
  #    We assume they are available in the namespace (no local definition).
  
  # row transform => sqrt(Mx) and col transform => sqrt(Ax)
  Xrw <- row_transform(Xp, build_sqrt_mult(Mx, rank_Mx, var_threshold, max_k, tol=tol))
  Xstar <- col_transform(Xrw, build_sqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol))
  
  Yrw <- row_transform(Yp, build_sqrt_mult(My, rank_My, var_threshold, max_k, tol=tol))
  Ystar <- col_transform(Yrw, build_sqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol))
  
  if(verbose) {
    message(sprintf("genplscorr: forming cross-product C, with p=%d x q=%d", p, q))
  }
  C <- crossprod(Xstar, Ystar)  # dimension p x q
  
  if(verbose) message("genplscorr: performing SVD on cross-product.")
  # For partial approach with big p,q, you might use RSpectra or PRIMME, but here base::svd is fine
  svdres <- tryCatch({
    if(ncomp >= min(p,q)) {
      base::svd(C, nu=min(p,q), nv=min(p,q))
    } else {
      tmp <- base::svd(C, nu=min(p,q), nv=min(p,q))
      tmp
    }
  }, error=function(e) {
    stop("genplscorr: SVD of cross-product C failed: ", e$message)
  })
  
  d_full <- svdres$d
  u_full <- svdres$u
  v_full <- svdres$v
  r0 <- min(length(d_full), ncomp)
  
  d_trunc <- d_full[seq_len(r0)]
  u_trunc <- u_full[, seq_len(r0), drop=FALSE]
  v_trunc <- v_full[, seq_len(r0), drop=FALSE]
  
  # factor scores in embedded domain (optionally):
  Tstar <- Xstar %*% u_trunc
  Ustar <- Ystar %*% v_trunc
  
  # map loadings back to original domain
  Ax_invsqrt <- build_invsqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
  Ay_invsqrt <- build_invsqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
  
  vx_emb <- u_trunc
  vy_emb <- v_trunc
  
  vx_orig <- Ax_invsqrt(vx_emb)
  vy_orig <- Ay_invsqrt(vy_emb)
  
  
  
  out <- cross_projector(
    vx = vx_orig,
    vy = vy_orig,
    preproc_x = px,
    preproc_y = py,
    classes   = c("genplscorr"),
    tilde_vx = vx_emb,
    tilde_vy = vy_emb,
    Tx = Tstar,
    Ty = Ustar,
    ncomp = r0,
    d_full = d_full,              # store full singular values
    rank_Mx = rank_Mx,
    rank_Ax = rank_Ax,
    rank_My = rank_My,
    rank_Ay = rank_Ay,
    var_threshold = var_threshold,
    max_k = max_k,
    tol = tol,
    verbose = verbose
  )
  out
}